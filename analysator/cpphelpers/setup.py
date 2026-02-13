from setuptools import setup, Extension
import numpy
setup(
    name="cpphelpers",
    version="1.0",
    ext_modules=[Extension("cpphelpers", sources=["cpphelpers.cpp"])],
    include_dirs=[numpy.get_include()]
)
