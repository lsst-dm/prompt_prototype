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

import collections.abc
import itertools
import logging
import os
import os.path
import tempfile

from lsst.utils import getPackageDir
from lsst.resources import ResourcePath
import lsst.afw.cameraGeom
from lsst.ctrl.mpexec import SimplePipelineExecutor
from lsst.daf.butler import Butler, CollectionType
import lsst.geom
from lsst.meas.algorithms.htmIndexer import HtmIndexer
import lsst.obs.base
import lsst.pipe.base

from .visit import Visit

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
    _COLLECTION_TEMPLATE = "templates"
    """The collection used for templates.
    """
    _COLLECTION_SKYMAP = "skymaps"
    """The collection used for skymaps.
    """

    def __init__(self, central_butler: Butler, image_bucket: str, instrument: str,
                 butler: Butler,
                 prefix: str = "gs://"):
        self.ip_apdb = os.environ["IP_APDB"]
        self.prefix = prefix
        self.central_butler = central_butler
        self.image_bucket = image_bucket
        # TODO: _download_store turns MWI into a tagged class; clean this up later
        if not prefix.startswith("file"):
            self._download_store = tempfile.TemporaryDirectory(prefix="holding-")
        else:
            self._download_store = None
        self.instrument = lsst.obs.base.Instrument.from_string(instrument)

        self.output_collection = self.instrument.makeCollectionName("prompt")

        self._init_local_butler(butler)
        self._init_ingester()

        define_visits_config = lsst.obs.base.DefineVisitsConfig()
        self.define_visits = lsst.obs.base.DefineVisitsTask(config=define_visits_config, butler=self.butler)

        # TODO DM-34098: note that we currently need to supply instrument here.
        # HACK: explicit collection gets around the fact that we don't have any
        # timestamp/exposure information in a form we can pass to the Butler.
        # This code will break once cameras start being versioned.
        self.camera = self.central_butler.get(
            "camera", instrument=self.instrument.getName(),
            collections=self.instrument.makeUnboundedCalibrationRunName()
        )
        self.skymap = self.central_butler.get("skyMap")

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
        self.instrument.register(butler.registry)

        # Refresh butler after configuring it, to ensure all required
        # dimensions and collections are available.
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

    def _predict_wcs(self, detector: lsst.afw.cameraGeom.Detector, visit: Visit) -> lsst.afw.geom.SkyWcs:
        """Calculate the expected detector WCS for an incoming observation.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            The detector for which to generate a WCS.
        visit : `Visit`
            Predicted observation metadata for the detector.

        Returns
        -------
        wcs : `lsst.afw.geom.SkyWcs`
            An approximate WCS for ``visit``.
        """
        boresight_center = lsst.geom.SpherePoint(visit.ra, visit.dec, lsst.geom.degrees)
        orientation = lsst.geom.Angle(visit.rot, lsst.geom.degrees)
        flip_x = True if self.instrument.getName() == "DECam" else False
        return lsst.obs.base.createInitialSkyWcsFromBoresight(boresight_center,
                                                              orientation,
                                                              detector,
                                                              flipX=flip_x)

    def _detector_bounding_circle(self, detector: lsst.afw.cameraGeom.Detector,
                                  wcs: lsst.afw.geom.SkyWcs
                                  ) -> (lsst.geom.SpherePoint, lsst.geom.Angle):
        # Could return a sphgeom.Circle, but that would require a lot of
        # sphgeom->geom conversions downstream. Even their Angles are different!
        """Compute a small sky circle that contains the detector.

        Parameters
        ----------
        detector : `lsst.afw.cameraGeom.Detector`
            The detector for which to compute an on-sky bounding circle.
        wcs : `lsst.afw.geom.SkyWcs`
            The conversion from detector to sky coordinates.

        Returns
        -------
        center : `lsst.geom.SpherePoint`
            The center of the bounding circle.
        radius : `lsst.geom.Angle`
            The opening angle of the bounding circle.
        """
        radii = []
        center = wcs.pixelToSky(detector.getCenter(lsst.afw.cameraGeom.PIXELS))
        for corner in detector.getCorners(lsst.afw.cameraGeom.PIXELS):
            radii.append(wcs.pixelToSky(corner).separation(center))
        return center, max(radii)

    def prep_butler(self, visit: Visit) -> None:
        """Prepare a temporary butler repo for processing the incoming data.

        Parameters
        ----------
        visit : `Visit`
            Group of snaps from one detector to prepare the butler for.
        """
        _log.info(f"Preparing Butler for visit '{visit}'")

        detector = self.camera[visit.detector]
        wcs = self._predict_wcs(detector, visit)
        center, radius = self._detector_bounding_circle(detector, wcs)

        with tempfile.NamedTemporaryFile(mode="w+b", suffix=".yaml") as export_file:
            with self.central_butler.export(filename=export_file.name, format="yaml") as export:
                self._export_refcats(export, center, radius)
                self._export_skymap_and_templates(export, center, detector, wcs)
                self._export_calibs(export, visit.detector, visit.filter)

                # CHAINED collections
                export.saveCollection(self.instrument.makeRefCatCollectionName())
                export.saveCollection(self._COLLECTION_TEMPLATE)
                export.saveCollection(self.instrument.makeUmbrellaCollectionName())

            self.butler.import_(filename=export_file.name,
                                directory=self.central_butler.datastore.root,
                                transfer="copy")

        self._prep_collections()
        self._prep_pipeline(visit)

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
        export.saveDatasets(self.central_butler.registry.queryDatasets(
            possible_refcats,
            collections=self.instrument.makeRefCatCollectionName(),
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
        # TODO: This exports the whole skymap, but we want to only export the
        # subset of the skymap that covers this data.
        # TODO: We only want to import the skymap dimension once in init,
        # otherwise we get a UNIQUE constraint error when prepping for the
        # second visit.
        export.saveDatasets(self.central_butler.registry.queryDatasets("skyMap",
                                                                       collections=self._COLLECTION_SKYMAP,
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
                                                                       collections=self._COLLECTION_TEMPLATE,
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
        export.saveDatasets(
            self.central_butler.registry.queryDatasets(
                ...,
                collections=self.instrument.makeCalibrationCollectionName(),
                where=calib_where),
            elements=[])  # elements=[] means do not export dimension records
        target_types = {CollectionType.CALIBRATION}
        for collection in self.central_butler.registry.queryCollections(...,
                                                                        collectionTypes=target_types):
            export.saveCollection(collection)

    def _prep_collections(self):
        """Pre-register output collections in advance of running the pipeline.
        """
        # NOTE: Because we receive a butler on init, we can't use this
        # prep_butler() because it takes a repo path.
        # butler = SimplePipelineExecutor.prep_butler(
        #     self.repo,
        #     inputs=[self.calibration_collection,
        #             self.instrument.makeDefaultRawIngestRunName(),
        #             'refcats'],
        #     output=self.output_collection)
        # The below is taken from SimplePipelineExecutor.prep_butler.
        # TODO DM-34202: save run collection in self for now, but we won't need
        # it when we no longer need to work around DM-34202.
        self.output_run = f"{self.output_collection}/{self.instrument.makeCollectionTimestamp()}"
        self.butler.registry.registerCollection(self.instrument.makeDefaultRawIngestRunName(),
                                                CollectionType.RUN)
        self.butler.registry.registerCollection(self.output_run, CollectionType.RUN)
        self.butler.registry.registerCollection(self.output_collection, CollectionType.CHAINED)
        collections = [self.instrument.makeUmbrellaCollectionName(),
                       self.instrument.makeDefaultRawIngestRunName(),
                       self.output_run]
        self.butler.registry.setCollectionChain(self.output_collection, collections)

        # Need to create a new butler with all the output collections.
        self.butler = Butler(butler=self.butler,
                             collections=[self.output_collection],
                             # TODO DM-34202: hack around a middleware bug.
                             run=None)

    def _prep_pipeline(self, visit: Visit) -> None:
        """Setup the pipeline to be run, based on the configured instrument and
        details of the incoming visit.

        Parameters
        ----------
        visit : Visit
            Group of snaps from one detector to prepare the pipeline for.

        Raises
        ------
        RuntimeError
            Raised if there is no AP pipeline file for this configuration.
            TODO: could be a good case for a custom exception here.
        """
        # TODO: We hacked the basepath in the Dockerfile so this works both in
        # development and in service container, but it would be better if there
        # were a path that's valid in both.
        ap_pipeline_file = os.path.join(getPackageDir("prompt_prototype"),
                                        "pipelines", self.instrument.getName(), "ApPipe.yaml")
        try:
            self.pipeline = lsst.pipe.base.Pipeline.fromFile(ap_pipeline_file)
        except FileNotFoundError:
            raise RuntimeError(f"No ApPipe.yaml defined for camera {self.instrument.getName()}")
        # TODO: Can we write to a configurable apdb schema, rather than
        # "postgres"?
        self.pipeline.addConfigOverride("diaPipe", "apdb.db_url",
                                        f"postgresql://postgres@{self.ip_apdb}/postgres")

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
        # TODO: this requires the service account to have the otherwise admin-ish
        # storage.buckets.get permission (DM-34188). Once that's resolved, see if
        # prompt-service can do without the "Storage Legacy Bucket Reader" role.
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

    def run_pipeline(self, visit: Visit, exposure_ids: set) -> None:
        """Process the received image(s).

        Parameters
        ----------
        visit : Visit
            Group of snaps from one detector to be processed.
        exposure_ids : `set`
            Identifiers of the exposures that were received.
        """
        # TODO: we want to define visits earlier, but we have to ingest a
        # faked raw file and appropriate SSO data during prep (and then
        # cleanup when ingesting the real data).
        try:
            self.define_visits.run({"instrument": self.instrument.getName(),
                                    "exposure": exp} for exp in exposure_ids)
        except lsst.daf.butler.registry.DataIdError as e:
            # TODO: a good place for a custom exception?
            raise RuntimeError("No data to process.") from e

        # TODO: can we move this from_pipeline call to prep_butler?
        where = f"detector={visit.detector} and exposure in ({','.join(str(x) for x in exposure_ids)})"
        executor = SimplePipelineExecutor.from_pipeline(self.pipeline, where=where, butler=self.butler)
        if len(executor.quantum_graph) == 0:
            # TODO: a good place for a custom exception?
            raise RuntimeError("No data to process.")
        # TODO DM-34202: hack around a middleware bug.
        executor.butler = Butler(butler=self.butler,
                                 collections=[self.output_collection],
                                 run=self.output_run)
        _log.info(f"Running '{self.pipeline._pipelineIR.description}' on {where}")
        # If this is a fresh (local) repo, then types like calexp,
        # *Diff_diaSrcTable, etc. have not been registered.
        result = executor.run(register_dataset_types=True)
        _log.info(f"Pipeline successfully run on {len(result)} quanta.")


def _query_missing_datasets(src_repo: Butler, dest_repo: Butler,
                            *args, **kwargs) -> collections.abc.Iterable[lsst.daf.butler.DatasetRef]:
    """Return datasets that are present in one repository but not another.

    Parameters
    ----------
    src_repo : `lsst.daf.butler.Butler`
        The repository in which a dataset must be present.
    dest_repo : `lsst.daf.butler.Butler`
        The repository in which a dataset must not be present.
    *args, **kwargs
        Parameters for describing the dataset query. They have the same
        meanings as the parameters of `lsst.daf.butler.Registry.queryDatasets`.

    Returns
    -------
    datasets : iterable [`lsst.daf.butler.DatasetRef`]
        The datasets that exist in ``src_repo`` but not ``dest_repo``.
    """
    try:
        # TODO: storing this set in memory may be a performance bottleneck.
        # In general, expect dest_repo to have more datasets than src_repo.
        known_datasets = set(dest_repo.registry.queryDatasets(*args, **kwargs))
    except lsst.daf.butler.registry.DataIdValueError as e:
        _log.debug("Pre-export query with args '%s, %s' failed with %s",
                   ", ".join(str(a) for a in args),
                   ", ".join(f"{k}={v}" for k, v in kwargs.items()),
                   e)
        # If dimensions are invalid, then *any* such datasets are missing.
        known_datasets = set()
    # If src_repo query fails, that invalidates this operation; raise it!
    return itertools.filterfalse(lambda ref: ref in known_datasets,
                                 src_repo.registry.queryDatasets(*args, **kwargs))
