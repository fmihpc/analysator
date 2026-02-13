from setuptools import setup, Extension
import numpy
setup(
    name="cpphelpers",
    version="1.0",
    description="Function to bind python and C++ together and speed pu things,requires numpy.",
    ext_modules=[Extension("cpphelpers", sources=["cpphelpers.cpp"])],
    include_dirs=[numpy.get_include()]
)
