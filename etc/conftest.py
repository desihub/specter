# This is a configuration file for pytest, and is set up to work with the
# pytest-astropy plugin.
#
# This *specific* file is intended for use with the Dockerfile in the same
# directory.
try:
    from pytest_astropy_header.display import PYTEST_HEADER_MODULES, TESTED_VERSIONS
    ASTROPY_HEADER = True
except ImportError:
    ASTROPY_HEADER = False


def pytest_configure(config):

    if ASTROPY_HEADER:
        config.option.astropy_header = True
        #
        # Customize the following lines to add/remove entries from the list of
        # packages for which version numbers are displayed when running the tests.
        #
        PYTEST_HEADER_MODULES.pop('Pandas', None)
        PYTEST_HEADER_MODULES.pop('h5py', None)
        PYTEST_HEADER_MODULES.pop('Matplotlib', None)
        PYTEST_HEADER_MODULES['Astropy'] = 'astropy'
        PYTEST_HEADER_MODULES['Numba'] = 'numba'
        PYTEST_HEADER_MODULES['fitsio'] = 'fitsio'
