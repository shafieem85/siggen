import pybind11
from distutils.core import setup, Extension


class get_pybind_include(object):
    """Helper class to determine the pybind11 include path
    The purpose of this class is to postpone importing pybind11
    until it is actually installed, so that the ``get_include()``
    method can be invoked. """

    def __str__(self):
        import pybind11
        return pybind11.get_include()
        
        

ext_modules = [
    Extension(
        'mjd_fieldgen',
        # Sort input source files to ensure bit-for-bit reproducible builds
        # (https://github.com/pybind/python_example/pull/53)
        sorted(['all.cpp']),
        include_dirs=[
            # Path to pybind11 headers
            get_pybind_include(),
        ],
        language='c++'
    ),
]

setup(
    name='mjd_fieldgen', # имя библиотеки собранной pybind11
    version='1.0.0',
    author='radforddc,DFed',
    author_email='random',
    description='pybind11 extension',
    ext_modules=ext_modules,
    requires=['pybind11']  # Указываем зависимость от pybind11
    # package_dir = {'': 'lib'}
)

