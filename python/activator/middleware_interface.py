# This file is part of prompt_prototype.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

__all__ = ["MiddlewareInterface"]

import logging
import os
import tempfile

from lsst.resources import ResourcePath
import lsst.afw.cameraGeom
from lsst.daf.butler import Butler, CollectionType
from lsst.meas.algorithms.htmIndexer import HtmIndexer
import lsst.obs.base
import lsst.geom

from .visit import Visit

PIPELINE_MAP = dict(
    BIAS="bias.yaml",
    DARK="dark.yaml",
    FLAT="flat.yaml",
)

_log = logging.getLogger("lsst." + __name__)
_log.setLevel(logging.DEBUG)


class MiddlewareInterface:
    """Interface layer between the Butler middleware and the prompt processing
    data handling system, to handle processing individual images.

    An instance of this class will accept an incoming group of single-detector
    snaps to process, using an instance-local butler repo. The instance can
    pre-load the necessary calibrations to process an incoming detector-visit,
    ingest the data when it is available, and run the difference imaging
    pipeline, all in that local butler.

    Parameters
    ----------
    central_butler : `lsst.daf.butler.Butler`
        Butler repo containing the calibration and other data needed for
        processing images as they are received. This butler must be created
        with the default instrument, skymap, and collections assigned.
    image_bucket : `str`
        Storage bucket where images will be written to as they arrive.
        See also ``prefix``.
    instrument : `str`
        Full class name of the instrument taking the data, for populating
        butler collections and dataIds. Example: "lsst.obs.lsst.LsstCam"
        TODO: this arg can probably be removed and replaced with internal
        use of the butler.
    butler : `lsst.daf.butler.Butler`
        Local butler to process data in and hold calibrations, etc.; must be
        writeable.
    prefix : `str`, optional
        URI identification prefix; prepended to ``image_bucket`` when
        constructing URIs to retrieve incoming files. The default is
        appropriate for use in the Google Cloud environment; typically only
        change this when running local tests.
    """
    def __init__(self, central_butler: Butler, image_bucket: str, instrument: str,
                 butler: Butler,
                 prefix: str = "gs://"):
        self.prefix = prefix
        self.central_butler = central_butler
        self.image_bucket = image_bucket
        # TODO: this member turns MWI into a tagged class; clean this up later
        if not prefix.startswith("file"):
            self._download_store = tempfile.TemporaryDirectory(prefix="holding-")
        else:
            self._download_store = None
        self.instrument = lsst.obs.base.utils.getInstrument(instrument)

        self._init_local_butler(butler)
        self._init_ingester()
        # TODO DM-34098: note that we currently need to supply instrument here.
        # HACK: explicit collection gets around the fact that we don't have any
        # timestamp/exposure information in a form we can pass to the Butler.
        # This code will break once cameras start being versioned.
        self.camera = self.central_butler.get(
            "camera", instrument=self.instrument.getName(),
            collections=self.instrument.makeCalibrationCollectionName("unbounded")
        )
        self.skymap = self.central_butler.get("skyMap")

        # self.r = self.src.registry
        self.calibration_collection = f"{instrument}/calib"

        # How much to pad the refcat region we will copy over.
        self.padding = 30*lsst.geom.arcseconds

    def _init_local_butler(self, butler: Butler):
        """Prepare the local butler to ingest into and process from.

        Parameters
        ----------
        butler : `lsst.daf.butler.Butler`
            Local butler to process data in and hold calibrations, etc.; must
            be writeable.
        """
        # TODO: Replace registring the instrument with importing a "base" repo
        # structure from an export.yaml file.
        # TODO: Cannot do this until we have a way to only extract calibs we
        # don't already have, otherwise we get a unique constraint error when
        # importing the export in prep_butler().
        # self.instrument.register(butler.registry)
        # Refresh butler after configuring it, to ensure all required
        # dimensions are available.
        butler.registry.refresh()
        self.butler = butler

    def _init_ingester(self):
        """Prepare the raw file ingester to receive images into this butler.
        """
        config = lsst.obs.base.RawIngestConfig()
        config.transfer = "copy"  # Copy files into the local butler.
        # TODO: Could we use the `on_ingest_failure` and `on_success` callbacks
        # to send information back to this interface?
        config.failFast = True  # We want failed ingests to fail immediately.
        self.rawIngestTask = lsst.obs.base.RawIngestTask(config=config,
                                                         butler=self.butler)

    def prep_butler(self, visit: Visit) -> None:
        """Prepare a temporary butler repo for processing the incoming data.

        Parameters
        ----------
        visit : Visit
            Group of snaps from one detector to prepare the butler for.
        """
        _log.info(f"Preparing Butler for visit '{visit}'")

        with tempfile.NamedTemporaryFile(mode="w+b", suffix=".yaml") as export_file:
            with self.central_butler.export(filename=export_file.name, format="yaml") as export:
                boresight_center = lsst.geom.SpherePoint(visit.ra, visit.dec, lsst.geom.degrees)
                orientation = lsst.geom.Angle(visit.rot, lsst.geom.degrees)
                detector = self.camera[visit.detector]
                # TODO: where do we get flipX from? See RFC-605
                wcs = lsst.obs.base.createInitialSkyWcsFromBoresight(boresight_center,
                                                                     orientation,
                                                                     detector,
                                                                     flipX=False)
                radii = []
                center = wcs.pixelToSky(detector.getCenter(lsst.afw.cameraGeom.PIXELS))
                for corner in detector.getCorners(lsst.afw.cameraGeom.PIXELS):
                    radii.append(wcs.pixelToSky(corner).separation(center))
                radius = max(radii)

                self._export_refcats(export, center, radius)
                self._export_skymap_and_templates(export, center, detector, wcs)
                self._export_calibs(export, visit.detector, visit.filter)

                # CHAINED collections
                export.saveCollection("refcats")

            self.butler.import_(filename=export_file.name,
                                directory=self.central_butler.datastore.root,
                                transfer="copy")

    def _export_refcats(self, export, center, radius):
        """Export the refcats for this visit from the central butler.

        Parameters
        ----------
        export : `Iterator[RepoExportContext]`
            Export context manager.
        center : `lsst.geom.SpherePoint`
            Center of the region to find refcat shards in.
        radius : `lst.geom.Angle`
            Radius to search for refcat shards in.
        """
        indexer = HtmIndexer(depth=7)
        shards = indexer.getShardIds(center, radius+self.padding)
        # getShardIds returns a tuple, the first item is the ids list.
        htm_where = f"htm7 in ({','.join(str(x) for x in shards[0])})"
        # Get shards from all refcats that overlap this detector.
        # TODO: `...` doesn't work for this queryDatasets call
        # currently, and we can't queryDatasetTypes in just the refcats
        # collection, so we have to specify a list here. Replace this
        # with another solution ASAP.
        possible_refcats = ["gaia", "panstarrs", "gaia_dr2_20200414", "ps1_pv3_3pi_20170110"]
        export.saveDatasets(self.central_butler.registry.queryDatasets(possible_refcats,
                                                                       collections="refcats",
                                                                       where=htm_where,
                                                                       findFirst=True))

    def _export_skymap_and_templates(self, export, center, detector, wcs):
        """Export the skymap and templates for this visit from the central
        butler.

        Parameters
        ----------
        export : `Iterator[RepoExportContext]`
            Export context manager.
        center : `lsst.geom.SpherePoint`
            Center of the region to load the skyamp tract/patches for.
        detector : `lsst.afw.cameraGeom.Detector`
            Detector we are loading data for.
        wcs : `lsst.afw.geom.SkyWcs`
            Rough WCS for the upcoming visit, to help finding patches.
        """
        export.saveDatasets(self.central_butler.registry.queryDatasets("skyMap",
                                                                       collections="skymaps",
                                                                       findFirst=True))
        # Getting only one tract should be safe: we're getting the
        # tract closest to this detector, so we should be well within
        # the tract bbox.
        tract = self.skymap.findTract(center)
        points = [center]
        for corner in detector.getCorners(lsst.afw.cameraGeom.PIXELS):
            points.append(wcs.pixelToSky(corner))
        patches = tract.findPatchList(points)
        patches_str = ','.join(str(p.sequential_index) for p in patches)
        template_where = f"patch in ({patches_str}) and tract={tract.tract_id}"
        # TODO: do we need to have the coadd name used in the pipeline
        # specified as a class kwarg, so that we only load one here?
        # TODO: alternately, we need to extract it from the pipeline? (best?)
        # TODO: alternately, can we just assume that there is exactly
        # one coadd type in the central butler?
        export.saveDatasets(self.central_butler.registry.queryDatasets("*Coadd",
                                                                       collections="templates",
                                                                       where=template_where))

    def _export_calibs(self, export, detector_id, filter):
        """Export the calibs for this visit from the central butler.

        Parameters
        ----------
        export : `Iterator[RepoExportContext]`
            Export context manager.
        detector_id : `int`
            Identifier of the detector to load calibs for.
        filter : `str`
            Physical filter name of the upcoming visit.
        """
        # TODO: we can't filter by validity range because it's not
        # supported in queryDatasets yet.
        calib_where = f"detector={detector_id} and physical_filter='{filter}'"
        export.saveDatasets(self.central_butler.registry.queryDatasets(
            ...,
            collections=self.instrument.makeCalibrationCollectionName(),
            where=calib_where))
        target_types = {CollectionType.CALIBRATION}
        for collection in self.central_butler.registry.queryCollections(...,
                                                                        collectionTypes=target_types):
            export.saveCollection(collection)

    def _download(self, remote):
        """Download an image located on a remote store.

        Parameters
        ----------
        remote : `lsst.resources.ResourcePath`
            The location from which to download the file. Must not be a
            file:// URI.

        Returns
        -------
        local : `lsst.resources.ResourcePath`
            The location to which the file has been downloaded.
        """
        local = ResourcePath(os.path.join(self._download_store.name, remote.basename()))
        local.transfer_from(remote, "copy")
        return local

    def ingest_image(self, oid: str) -> None:
        """Ingest an image into the temporary butler.

        Parameters
        ----------
        oid : `str`
            Google storage identifier for incoming image, relative to the
            image bucket.
        """
        _log.info(f"Ingesting image id '{oid}'")
        file = ResourcePath(f"{self.prefix}{self.image_bucket}/{oid}")
        if not file.isLocal:
            # TODO: RawIngestTask doesn't currently support remote files.
            file = self._download(file)
        result = self.rawIngestTask.run([file])
        # We only ingest one image at a time.
        # TODO: replace this assert with a custom exception, once we've decided
        # how we plan to handle exceptions in this code.
        assert len(result) == 1, "Should have ingested exactly one image."
        _log.info("Ingested one %s with dataId=%s", result[0].datasetType.name, result[0].dataId)

    def run_pipeline(self, visit: Visit, snaps: set) -> None:
        """Process the received image.

        Parameters
        ----------
        visit : Visit
            Group of snaps from one detector to be processed.
        snaps : `set`
            Identifiers of the snaps that were received.
            TODO: I believe this is unnecessary because it should be encoded
            in the `visit` object, but we'll have to test how that works once
            we implemented this with actual data.
        """
        pipeline = PIPELINE_MAP[visit.kind]
        _log.info(f"Running pipeline {pipeline} on visit '{visit}', snaps {snaps}")
        cmd = [
            "echo",
            "pipetask",
            "run",
            "-b",
            self.butler,
            "-p",
            pipeline,
            "-i",
            f"{self.instrument}/raw/all",
        ]
        _log.debug("pipetask command line: %s", cmd)
        # subprocess.run(cmd, check=True)


def filter_calibs(dataset_ref, visit_info):
    for dimension in ("instrument", "detector", "physical_filter"):
        if dimension in dataset_ref.dataId:
            if dataset_ref.dataId[dimension] != visit_info[dimension]:
                return False
    return True
