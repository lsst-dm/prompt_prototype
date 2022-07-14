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

import itertools
import tempfile
import os.path
import unittest
import unittest.mock

import astropy.coordinates
import astropy.units as u

import astro_metadata_translator
from lsst.daf.butler import Butler, DataCoordinate
from lsst.obs.base.formatters.fitsExposure import FitsImageFormatter
from lsst.obs.base.ingest import RawFileDatasetInfo, RawFileData
import lsst.resources

from activator.visit import Visit
from activator.middleware_interface import MiddlewareInterface, _query_missing_datasets

# The short name of the instrument used in the test repo.
instname = "DECam"
# Full name of the physical filter for the test file.
filter = "g DECam SDSS c0001 4720.0 1520.0"


def fake_file_data(filename, dimensions, instrument, visit):
    """Return file data for a mock file to be ingested.

    Parameters
    ----------
    filename : `str`
        Full path to the file to mock. Can be a non-existant file.
    dimensions : `lsst.daf.butler.DimensionsUniverse`
        The full set of dimensions for this butler.
    instrument : `lsst.obs.base.Instrument`
        The instrument the file is supposed to be from.
    visit : `Visit`
        Group of snaps from one detector to be processed.

    Returns
    -------
    data_id, file_data, : `DataCoordinate`, `RawFileData`
        The id and descriptor for the mock file.
    """
    data_id = DataCoordinate.standardize({"exposure": 1,
                                          "detector": visit.detector,
                                          "instrument": instrument.getName()},
                                         universe=dimensions)

    time = astropy.time.Time("2015-02-18T05:28:18.716517500", scale="tai")
    obs_info = astro_metadata_translator.makeObservationInfo(
        instrument=instrument.getName(),
        datetime_begin=time,
        datetime_end=time + 30*u.second,
        exposure_id=1,
        visit_id=1,
        boresight_rotation_angle=astropy.coordinates.Angle(visit.rot*u.degree),
        boresight_rotation_coord='sky',
        tracking_radec=astropy.coordinates.SkyCoord(visit.ra, visit.dec, frame="icrs", unit="deg"),
        observation_id="1",
        physical_filter=filter,
        exposure_time=30.0*u.second,
        observation_type="science")
    dataset_info = RawFileDatasetInfo(data_id, obs_info)
    file_data = RawFileData([dataset_info],
                            lsst.resources.ResourcePath(filename),
                            FitsImageFormatter,
                            instrument)
    return data_id, file_data


class MiddlewareInterfaceTest(unittest.TestCase):
    """Test the MiddlewareInterface class with faked data.
    """
    @classmethod
    def setUpClass(cls):
        cls.env_patcher = unittest.mock.patch.dict(os.environ,
                                                   {"IP_APDB": "localhost"})
        cls.env_patcher.start()

        super().setUpClass()

    @classmethod
    def tearDownClass(cls):
        super().tearDownClass()

        cls.env_patcher.stop()

    def setUp(self):
        data_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)), "data")
        central_repo = os.path.join(data_dir, "central_repo")
        central_butler = Butler(central_repo,
                                instrument=instname,
                                skymap="deepCoadd_skyMap",
                                collections=[f"{instname}/defaults"],
                                writeable=False)
        instrument = "lsst.obs.decam.DarkEnergyCamera"
        self.input_data = os.path.join(data_dir, "input_data")
        # Have to preserve the tempdir, so that it doesn't get cleaned up.
        self.repo = tempfile.TemporaryDirectory()
        self.butler = Butler(Butler.makeRepo(self.repo.name), writeable=True)
        self.interface = MiddlewareInterface(central_butler, self.input_data, instrument, self.butler,
                                             prefix="file://")

        # coordinates from DECam data in ap_verify_ci_hits2015 for visit 411371
        ra = 155.4702849608958
        dec = -4.950050405424033
        # DECam has no rotator; instrument angle is 90 degrees in our system.
        rot = 90.
        self.next_visit = Visit(instrument,
                                detector=56,
                                group=1,
                                snaps=1,
                                filter=filter,
                                ra=ra,
                                dec=dec,
                                rot=rot,
                                kind="SURVEY")
        self.logger_name = "lsst.activator.middleware_interface"

    def tearDown(self):
        super().tearDown()
        # TemporaryDirectory warns on leaks
        self.repo.cleanup()

    def test_init(self):
        """Basic tests of the initialized interface object.
        """
        # Ideas for things to test:
        # * On init, does the right kind of butler get created, with the right
        #   collections, etc?
        # * On init, is the local butler repo purely in memory?

        # Check that the butler instance is properly configured.
        instruments = list(self.butler.registry.queryDimensionRecords("instrument"))
        self.assertEqual(instname, instruments[0].name)

        # Check that the ingester is properly configured.
        self.assertEqual(self.interface.rawIngestTask.config.failFast, True)
        self.assertEqual(self.interface.rawIngestTask.config.transfer, "copy")

    def _check_imports(self, butler, detector):
        """Test that the butler has the expected supporting data.
        """
        self.assertEqual(butler.get('camera',
                                    instrument=instname,
                                    collections=[f"{instname}/calib/unbounded"]).getName(), instname)

        # Check that the right skymap is in the chained output collection.
        self.assertTrue(
            butler.datasetExists("skyMap",
                                 # TODO: we shouldn't need skymap here?
                                 skymap="deepCoadd_skyMap",
                                 collections=self.interface.output_collection)
        )

        # check that we got appropriate refcat shards
        loaded_shards = list(butler.registry.queryDataIds("htm7",
                                                          datasets="gaia",
                                                          collections="refcats"))

        # These shards were identified by plotting the objects in each shard
        # on-sky and overplotting the detector corners.
        # TODO DM-34112: check these shards again with some plots, once I've
        # determined whether ci_hits2015 actually has enough shards.
        expected_shards = [157394, 157401, 157405]
        self.assertEqual(expected_shards, [x['htm7'] for x in loaded_shards])
        # Check that the right calibs are in the chained output collection.
        try:
            self.assertTrue(
                butler.datasetExists('cpBias', detector=detector, instrument='DECam',
                                     collections="DECam/calib/20150218T000000Z")
                # TODO: Have to use the exact run collection, because we can't
                # query by validity range.
                # collections=self.interface.output_collection)
            )
        except LookupError:
            self.fail("Bias file missing from local butler.")
        try:
            self.assertTrue(
                butler.datasetExists('cpFlat', detector=detector, instrument='DECam',
                                     physical_filter=filter,
                                     collections="DECam/calib/20150218T000000Z")
                # TODO: Have to use the exact run collection, because we can't
                # query by validity range.
                # collections=self.interface.output_collection)
            )
        except LookupError:
            self.fail("Flat file missing from local butler.")

        # Check that the right templates are in the chained output collection.
        # Need to refresh the butler to get all the dimensions/collections.
        butler.registry.refresh()
        for patch in (464, 465, 509, 510):
            butler.datasetExists('deepCoadd', tract=0, patch=patch, band="g",
                                 # TODO: we shouldn't need skymap here?
                                 skymap="deepCoadd_skyMap",
                                 collections=self.interface.output_collection)

    def test_prep_butler(self):
        """Test that the butler has all necessary data for the next visit.
        """
        self.interface.prep_butler(self.next_visit)

        self._check_imports(self.butler, detector=56)

        # Check that we configured the right pipeline.
        self.assertEqual(self.interface.pipeline._pipelineIR.description,
                         "End to end Alert Production pipeline specialized for HiTS-2015")

    def test_prep_butler_twice(self):
        """prep_butler should have the correct calibs (and not raise an
        exception!) on a second run with the same, or a different detector.
        This explicitly tests the "you can't import something that's already
        in the local butler" problem that's related to the "can't register
        the skymap in init" problem.
        """
        self.interface.prep_butler(self.next_visit)
        # TODO: update next_visit with a new group number
        self.interface.prep_butler(self.next_visit)
        self._check_imports(self.butler, detector=56)

    def test_ingest_image(self):
        filename = "fakeRawImage.fits"
        filepath = os.path.join(self.input_data, filename)
        data_id, file_data = fake_file_data(filepath,
                                            self.butler.dimensions,
                                            self.interface.instrument,
                                            self.next_visit)
        with unittest.mock.patch.object(self.interface.rawIngestTask, "extractMetadata") as mock:
            mock.return_value = file_data
            self.interface.ingest_image(filename)

            datasets = list(self.butler.registry.queryDatasets('raw',
                                                               collections=[f'{instname}/raw/all']))
            self.assertEqual(datasets[0].dataId, data_id)
            # TODO: After raw ingest, we can define exposure dimension records
            # and check that the visits are defined

    def test_ingest_image_fails_missing_file(self):
        """Trying to ingest a non-existent file should raise.

        NOTE: this is currently a bit of a placeholder: I suspect we'll want to
        change how errors are handled in the interface layer, raising custom
        exceptions so that the activator can deal with them better. So even
        though all this is demonstrating is that if the file doesn't exist,
        rawIngestTask.run raises FileNotFoundError and that gets passed up
        through ingest_image(), we'll want to have a test of "missing file
        ingestion", and this can serve as a starting point.
        """
        filename = "nonexistentImage.fits"
        filepath = os.path.join(self.input_data, filename)
        data_id, file_data = fake_file_data(filepath,
                                            self.butler.dimensions,
                                            self.interface.instrument,
                                            self.next_visit)
        with unittest.mock.patch.object(self.interface.rawIngestTask, "extractMetadata") as mock, \
                self.assertRaisesRegex(FileNotFoundError, "Resource at .* does not exist"):
            mock.return_value = file_data
            self.interface.ingest_image(filename)
        # There should not be any raw files in the registry.
        datasets = list(self.butler.registry.queryDatasets('raw',
                                                           collections=[f'{instname}/raw/all']))
        self.assertEqual(datasets, [])

    def test_run_pipeline(self):
        """Test that running the pipeline uses the correct arguments.
        We can't run an actual pipeline because raw/calib/refcat/template data
        are all zeroed out.
        """
        # Have to setup the data so that we can create the pipeline executor.
        self.interface.prep_butler(self.next_visit)
        filename = "fakeRawImage.fits"
        filepath = os.path.join(self.input_data, filename)
        data_id, file_data = fake_file_data(filepath,
                                            self.butler.dimensions,
                                            self.interface.instrument,
                                            self.next_visit)
        with unittest.mock.patch.object(self.interface.rawIngestTask, "extractMetadata") as mock:
            mock.return_value = file_data
            self.interface.ingest_image(filename)

        with unittest.mock.patch("activator.middleware_interface.SimplePipelineExecutor.run") as mock_run:
            self.interface.run_pipeline(self.next_visit, {1})
        mock_run.assert_called_once_with(register_dataset_types=True)

    def test_run_pipeline_empty_quantum_graph(self):
        """Test that running a pipeline that results in an empty quantum graph
        (because the exposure ids are wrong), raises.
        """
        # Have to setup the data so that we can create the pipeline executor.
        self.interface.prep_butler(self.next_visit)
        filename = "fakeRawImage.fits"
        filepath = os.path.join(self.input_data, filename)
        data_id, file_data = fake_file_data(filepath,
                                            self.butler.dimensions,
                                            self.interface.instrument,
                                            self.next_visit)
        with unittest.mock.patch.object(self.interface.rawIngestTask, "extractMetadata") as mock:
            mock.return_value = file_data
            self.interface.ingest_image(filename)

        with self.assertRaisesRegex(RuntimeError, "No data to process"):
            self.interface.run_pipeline(self.next_visit, {2})

    def test_query_missing_datasets(self):
        """Test that query_missing_datasets provides the correct values.
        """
        # Much easier to create DatasetRefs with a real repo.
        butler = self.interface.central_butler
        dtype = butler.registry.getDatasetType("cpBias")
        data1 = lsst.daf.butler.DatasetRef(dtype, {"instrument": "DECam", "detector": 5})
        data2 = lsst.daf.butler.DatasetRef(dtype, {"instrument": "DECam", "detector": 25})
        data3 = lsst.daf.butler.DatasetRef(dtype, {"instrument": "DECam", "detector": 42})

        for src, existing in itertools.product([set(), {data1, data2}, {data1, data2, data3}], repeat=2):
            diff = src - existing
            src_butler = unittest.mock.Mock(**{"registry.queryDatasets.return_value": src})
            existing_butler = unittest.mock.Mock(**{"registry.queryDatasets.return_value": existing})

            with self.subTest(src=sorted(ref.dataId["detector"] for ref in src),
                              existing=sorted(ref.dataId["detector"] for ref in existing)):
                result = set(_query_missing_datasets(src_butler, existing_butler,
                                                     "cpBias", instrument="DECam"))
                src_butler.registry.queryDatasets.assert_called_once_with("cpBias", instrument="DECam")
                existing_butler.registry.queryDatasets.assert_called_once_with("cpBias", instrument="DECam")
                self.assertEqual(result, diff)

    def test_query_missing_datasets_nodim(self):
        """Test that query_missing_datasets provides the correct values when
        the destination repository is missing not only datasets, but the
        dimensions to define them.
        """
        # Much easier to create DatasetRefs with a real repo.
        butler = self.interface.central_butler
        dtype = butler.registry.getDatasetType("skyMap")
        data1 = lsst.daf.butler.DatasetRef(dtype, {"skymap": "mymap"})

        src_butler = unittest.mock.Mock(**{"registry.queryDatasets.return_value": {data1}})
        existing_butler = unittest.mock.Mock(
            **{"registry.queryDatasets.side_effect":
               lsst.daf.butler.registry.DataIdValueError(
                   "Unknown values specified for governor dimension skymap: {'mymap'}")
               })

        result = set(_query_missing_datasets(src_butler, existing_butler, "skyMap", ..., skymap="mymap"))
        src_butler.registry.queryDatasets.assert_called_once_with("skyMap", ..., skymap="mymap")
        self.assertEqual(result, {data1})
