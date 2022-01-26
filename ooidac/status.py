import os
import json
import uuid
from datetime import datetime as dt
from hashlib import md5

class Status(object):
    """A status keeping object that writes to JSON for caching"""
    def __init__(self, status_path):
        """Writes out the status.json file"""
        self.path = status_path
        if not os.path.exists(self.path):
            now_tstr = dt.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
            self.info = {
                "trajectory_name": "",
                "history": "{:s}: dataset created.".format(now_tstr),
                "date_created": now_tstr,
                "date_modified": now_tstr,
                "date_issued": now_tstr,
                "version": "1.0",
                "uuid": "",
                #"raw_directory": raw_data_path,
                #"nc_directory": nc_path,
                "next_profile_id": None, "files_processed": [],
                "profiles_created": [], "profiles_uploaded": [],
                "profile_to_data_map": []
            }
        else:
            with open(self.path, 'r') as fid:
                self.info = json.load(fid)

    def add_nc(self, nc_file, src_file, overwrite=False):
        """adds created profile nc files to status"""
        #nc_fn = os.path.basename(nc_file)
        nc_fn = nc_file
        data_map = [nc_file, src_file]
        file_exists = os.path.exists(nc_file)
        file_included = nc_fn in self.info['profiles_created']
        add_file = file_exists and (not(file_included) or overwrite)
        if overwrite and file_included:
            self.info['profiles_created'].remove(nc_fn)
            self.info['profile_to_data_map'].remove(data_map)
        if add_file:
            self.info['profiles_created'].append(nc_fn)
            self.info['profile_to_data_map'].append(data_map)
            self.write()

    def add_src(self, src_file, overwrite=False):
        """adds processed source data files to status"""
        #src_fn = os.path.basename(src_file)
        src_fn = src_file
        file_exists = os.path.exists(src_file)
        file_included = src_fn in self.info['files_processed']
        add_file = file_exists and (not (file_included) or overwrite)
        if overwrite and file_included:
            self.info['files_processed'].remove(src_fn)
        if add_file:
            self.info['files_processed'].append(src_fn)
            self.write()

    def add_upload(self, nc_file, overwrite=False):
        """adds uploaded nc data files to status"""
        #nc_fn = os.path.basename(nc_file)
        nc_fn = nc_file
        file_included = nc_fn in self.info['profiles_uploaded']
        add_file = not (file_included) or overwrite
        if overwrite and file_included:
            self.info['profiles_uploaded'].remove(nc_fn)
        if add_file:
            self.info['profiles_uploaded'].append(nc_fn)
            self.write()

    def update_history(self, addtl_msg):
        """updates the history message in status.json"""
        now_tstr = dt.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        self.info['history'] += "\n{:s} {:s}".format(nowtstr, addtl_msg)
        self.info['date_modified'] = now_tstr
        self.write()

    @property
    def processed(self):
        """a list of the basenames of the processed source files"""
        return list(map(os.path.basename, self.info['files_processed']))

    @property
    def created(self):
        """a list of the basenames of the profile nc files created"""
        return list(map(os.path.basename, self.info['profiles_created']))

    @property
    def uploaded(self):
        """a list of the basenames of the profile nc files uploaded"""
        return list(map(os.path.basename, self.info['profiles_uploaded']))

    @property
    def data_map(self):
        """a list of the mapped basenames of profile nc files to source data files"""
        return list(map(self._twocolbasename, self.info['profile_to_data_map']))

    def _twocolbasename(self, twocols):
        """a special function to map to the data_map property"""
        return [os.path.basename(twocols[0]), os.path.basename(twocols[1])]

    @property
    def next_profile_id(self):
        """the last used profile ID"""
        return self.info['next_profile_id']

    def update_major_version(self):
        """steps major version to a higher number and updates the uuid"""
        maj_ver = int(self.info['version'].split('.')[0])
        maj_ver += 1
        self.version = "{:d}.0".format(maj_ver)
        self.write()

    def update_minor_version(self):
        """steps minor version to a higher number and updates the uuid"""
        maj_ver = self.info['version'].split('.')[0]
        min_ver = int(self.info['version'].split('.')[1])
        min_ver += 1
        self.version = "{:s}.{:d}".format(maj_ver, min_ver)
        self.write()

    def update_modified_date(self):
        """updates the modified date to the current date and time"""
        now_tstr = dt.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
        self.info['date_modified'] = now_tstr
        self.write()

    @property
    def version(self):
        return self.info['version']

    @version.setter
    def version(self, value):
        self.info['version'] = value
        self.info['uuid'] = create_uuid(self.trajectory, value)
        self.write()

    @property
    def trajectory(self):
        """returns the trajectory ID"""
        return self.info['trajectory_name']

    @trajectory.setter
    def trajectory(self, value):
        """sets the trajectory ID"""
        self.info['trajectory_name'] = value
        self.info['uuid'] = create_uuid(value, self.info['version'])
        self.write()

    def write(self):
        """write out the json file status"""
        with open(self.path, 'w') as fid:
            json.dump(self.info, fid, indent=2)


def create_uuid(trajectory_name, version):
    """creates a unique ID with the uuid library
    from the trajectory name and run version"""
    hashable_str = trajectory_name + ";" + str(version)
    md5hash = md5(bytes(hashable_str, "utf-8"))
    return str(uuid.UUID(hex=md5hash.hexdigest()))