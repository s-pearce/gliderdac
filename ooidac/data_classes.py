import os
import logging
import numpy as np
from ooidac.readers.slocum import parse_dba
from ooidac.utilities import cluster_index

logger = logging.getLogger(os.path.basename(__name__))


class GliderDataParticle(object):
    def __init__(self):
        pass


class GliderData(object):
    """

    """
    def __init__(self, metadata, sensor_names, sensors, data):
        self.file_metadata = metadata
        # for key in metadata:
        #     self.__setattr__(key, metadata[key])
        if 'source_file' in metadata:
            self.source_file = metadata['source_file']
        else:
            self.source_file = ""
        self._data = data
        self._sensor_names = sensor_names
        self.sensors = sensors
        self.N = len(self._data)
        # config is meant to be a container for attaching any configuration
        # object data required to keep with the glider data
        # E.g. dictionary of deployment info, list of special variables, etc.
        self.config = None
        self.timesensorname = None
        self.depthsensorname = None
        self.pressuresensorname = None
        self.ballastname = None
        self.scitimesensorname = None
        self.depth = None
        self.ts = None
        if self.N > 0:
            self.set_ts()
            self.set_depth()

    # ToDo: fix this function or go back to simple call to _get_dataparticle
    def __getitem__(self, item):
        return_item = None
        if isinstance(item, str):
            return_item = self._get_dataparticle(item)
        elif len(item) == 2 and np.array(item).dtype == np.dtype('O'):
            try:
                return_item = self.slicedata(sensors=item[0], indices=item[1])
            except:
                logger.warning("Slicing index of this type not permitted")
        elif len(item) > 2 and isinstance(item, (list, np.ndarray)):
            dtype = np.array(item).dtype
            if dtype.str.startswith('<U'):    # case of multiple
                return_item = self.slicedata(sensors=item)
            elif dtype == np.dtype('int32'):  # case of multiple integers
                return_item = self.slicedata(indices=item)
        return return_item

    def __setitem__(self, key, sensor_particle):
        # if this data particle is already there (with the same name) don't
        # continue
        self.add_data(sensor_particle, key)

    def __len__(self):
        """Returns a value when the len function is used"""
        return self.N

    def __repr__(self):
        """Return value when the repr function is called"""
        return "<GliderData({:d} x {:d})>".format(self.m, self.N)

    @property
    def m(self):
        return len(self.sensor_names)

    @property
    def sensor_names(self):
        # This is a wrapper to the private attribute just so a new sensor
        # name cannot be added without adding the corresponding data to the end
        # of the self._data array.
        return self._sensor_names

    def add_data(self, sensor_particle, key=None):
        if not isinstance(sensor_particle, dict):
            logger.warning(
                "added data should be a dictionary with keys 'attrs', 'data', "
                "and 'sensor_name'.")
            return

        # check that the sensor particle meets the format.
        keys = list(sensor_particle.keys())
        keys.sort()
        for attr in ['attrs', 'data', 'sensor_name']:
            if attr not in keys:
                logger.warning('sensor_particle must have the attributes '
                               '"attrs", "data", and "sensor_name". Data not '
                               'added.')
                return

        # if this data particle is already there (with the same name) don't
        # continue
        if key and key != sensor_particle['sensor_name']:
            logger.warning(
                'New data {:s} does not match sensor_name {:s}'.format(
                    key, sensor_particle['sensor_name'])
            )
            return
        elif not key:
            key = sensor_particle['sensor_name']

        if key in self.sensor_names:
            logger.warning((
                'Data already exists, Not adding new data {:s}.').format(key)
            )
            return

        if len(sensor_particle['data']) != self.N:
            logger.warning(
                ('Data in added sensor {:s} is not the same '
                 'length as the rest of the data').format(key)
            )
            return
        data = sensor_particle.pop('data')
        self._data = np.append(self._data, data.reshape((self.N, 1)), axis=1)
        self.sensors[key] = sensor_particle
        self._sensor_names.append(key)

    def getdata(self, item):
        if item in self._sensor_names:
            idx = self._sensor_names.index(item)
            return self._data[:, idx]
        else:
            raise SensorError("Sensor {:s} is not available".format(item))

    def getdataslice(self, items):
        if isinstance(items, str):
            items = [items]
        idxs = []
        for item in items:
            if item in self._sensor_names:
                idx = self._sensor_names.index(item)
                idxs.append(idx)
            else:
                raise SensorError("Sensor {:s} is not available".format(item))
        if len(idxs) == 1:
            idxs = idxs[0]
        return self._data[:, idxs]

    def update_data(self, items, row_indices, values):
        row_indices = np.atleast_1d(row_indices)
        col_idxs = []
        for item in items:
            if item in self._sensor_names:
                idx = self._sensor_names.index(item)
                col_idxs.append(idx)
            else:
                raise SensorError("Sensor {:s} is not available".format(item))
        col_idxs = np.atleast_1d(col_idxs)
        # Don't want a try statement here, I want the np.array error to raise
        # if `values` does not fit into the indices given
        self._data[row_indices.reshape(len(row_indices), 1), col_idxs] = values

    def _get_dataparticle(self, item):
        if item in self._sensor_names:
            data_particle = self.sensors[item].copy()
            idx = self._sensor_names.index(item)
            data_particle['data'] = self._data[:, idx]
            return_item = data_particle
        else:
            # return_item = None
            raise SensorError("Sensor {:s} is not available".format(item))
        return return_item

    def slicedata(self, sensors=None, indices=None):
        if sensors is None and indices is None:
            return

        if sensors is not None:
            sensor_names = sensors
            sensor_defs = {}
            col_inds = []
            for sensor in sensors:
                sensor_defs[sensor] = self.sensors[sensor].copy()
                col_inds.append(self._sensor_names.index(sensor))
        else:
            sensor_names = self._sensor_names.copy()
            sensor_defs = self.sensors.copy()
            col_inds = slice(None)

        if indices is not None and sensors is not None:
            row_inds = np.array(indices)[:, np.newaxis]
        elif indices is not None and sensors is None:
            row_inds = np.array(indices)
        else:
            row_inds = slice(None)

        data = self._data[row_inds, col_inds]
        depth = self.depth[row_inds]
        new_gdata_slice = GliderData(
            self.file_metadata, sensor_names, sensor_defs, data)
        new_gdata_slice.set_depth(depth=depth)
        return new_gdata_slice

    def set_ts(self, timesensor=None):
        """Set the .ts timestamp attribute with the desired time sensor.
        Default tries `m_present_time`

        :param timesensor: Optional, name of time sensor to use
        """
        if timesensor is None:
            timesensor = 'm_present_time'
        try:
            self.ts = self.getdata(timesensor)
        except SensorError:
            logging.warning('Time sensor {:s} not found. {:s}'.format(
                timesensor, self.source_file
            ))
            self.ts = None
            return

    def set_depth(self, depthsensor=None, interpolate=True, depth=None):
        """Set the .depth attribute with the desired depth sensor.
        Default tries `m_depth`

        In case the GliderData instance is a slice from a larger GliderData
        instance, there is a chance that the slice selected has too few data
        points from the larger instance.  So in the case of a slice,
        this method has the parameter `depth` where a slice of a larger
        interpolated depth can be passed.

        :param depthsensor: Optional Str, name of depth sensor to use.
            default is `m_depth`
        :param interpolate: Optional Bool, interpolate the missing depth values
            default is True.
        :param depth: Optional Numpy array, array of depth values to use
            default is None
        """
        if depth is not None:
            self.depth = depth
            return
        if depthsensor is None:
            depthsensor = 'm_depth'
        try:
            depth = self.getdata(depthsensor)
        except SensorError:
            logging.warning('Depth sensor {:s} not found. {:s}'.format(
                depthsensor, self.source_file
            ))
            self.depth = None
            return
        if (
                interpolate and np.any(np.isfinite(depth))
                and self.ts is not None
        ):
            finites = np.flatnonzero(np.isfinite(depth))
            self.depth = np.interp(
                self.ts, self.ts[finites], depth[finites],
                left=depth[finites[0]], right=depth[finites[-1]])
        else:
            self.depth = depth


class DbaData(GliderData):
    """

    """
    def __init__(self, dba_file):
        # self.file_metadata = None
        # self._data = np.array([])
        dba = parse_dba(dba_file)
        if dba is None:
            return
        else:
            GliderData.__init__(
                self,
                dba['header'], dba['sensor_names'],
                dba['sensor_defs'], dba['data']
            )
        self.underwater_indices = None
        self.pre_dive_indices = None
        self.post_dive_indices = None
        self.surface_indices = None
        if self.N > 0:
            self.get_indices()

    def get_indices(self):
        """ Determine the indices for when the glider is underwater, at the
        surface, and if the glider dives, the pre-dive and post dive indices.
        """
        if self.depth is None:
            return
        underwater_indices = np.flatnonzero(self.depth > 2.0)
        # max_depth = self.depth.max()
        clusters, cluster_ids = cluster_index(underwater_indices, ids=True)
        for cluster_id in cluster_ids:
            cluster = np.flatnonzero(clusters == cluster_id)
            if len(cluster) < 5:
                clusters = np.delete(clusters, cluster)
                underwater_indices = np.delete(underwater_indices, cluster)
        self.underwater_indices = underwater_indices
        self.surface_indices = np.setdiff1d(
            np.arange(self.N), underwater_indices)  # order matters
        if len(underwater_indices > 0):
            self.pre_dive_indices = np.arange(underwater_indices[0])
            self.post_dive_indices = np.arange(
                underwater_indices[-1]+1, self.N)


class SensorError(KeyError):
    pass
