#!/usr/bin/python
#
# Fourier transform a scalar variable, or one component of a vector, in a vlsv file
# (entiere space)
#
import pytools as pt
import numpy
import pylab
import sys

if len(sys.argv) < 3:
    sys.stderr.write("Syntax: fourier_transform_field.py <vlsvfile> <field> [<component>]\n")
    sys.exit()

filename = sys.argv[1]
variable = sys.argv[2]

# Open file
f = pt.vlsvfile.VlsvReader(filename)

# Read our field and cellIDs to sort it
cellids = f.read_variable("CellID")
a = f.read_variable(variable)

if len(a.shape) > 1:
    if len(sys.argv) > 3:
        a = a[:,int(sys.argv[3])]
    else:
        a = numpy.sqrt( numpy.sum(a*a,1))

xsize = f.read_parameter("xcells_ini")
ysize = f.read_parameter("ycells_ini")

# Sort the Field to be a proper 2D numpy-array
a = a[cellids.argsort()].reshape([ysize,xsize])

# 2D fourier-transform it.
af = numpy.fft.fftshift(numpy.fft.fft2(a))


# Plot it!
pylab.plt.subplot().pcolormesh(numpy.log(abs(af)))
pylab.plt.show()
