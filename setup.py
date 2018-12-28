import setuptools

def configure():
# Initialize the setup kwargs
    kwargs = {
            'name': 'h5seis',
            'version': '0.0a0',
            'author': 'Malcolm White',
            'author_email': 'malcolm.white@usc.edu',
            'maintainer': 'Malcolm White',
            'maintainer_email': 'malcolm.white@usc.edu',
            'url': 'http://malcolmw.github.io/h5seis',
            'description': 'HDF5-based seismic data storage.',
            'download_url': 'https://github.com/malcolmw/h5seis.git',
            'platforms': ['linux'],
            'requires': ['numpy', 'obspy', 'pandas', 'scipy'],
            'packages': ['h5seis']
            }
    return(kwargs)

if __name__ == '__main__':
    kwargs = configure()
    setuptools.setup(**kwargs)
