from setuptools import setup, find_packages
from ooidac import __version__


def readme():
    with open('README.md') as f:
        return f.read()


reqs = [line.strip() for line in open('requirements.txt')]

setup(
    name="gliderdac",
    version=__version__,
    description="Python config_file_conversions for writing profile NetCDF "
                "files from raw Slocum glider data files for upload to the "
                "National Glider Data Assembly Center",
    long_description=readme(),
    packages=find_packages(),
    author="Stuart Pearce, John Kerfoot",
    author_email="stuart.pearce@oregonstate.edu",
    license="GPLv3",
    install_requires=reqs,
    tests_require=['pytest'],
)
