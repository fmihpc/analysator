import pytools as pt

f = pt.vlsvfile.VlasiatorReader('fullf.0000002.vlsv')
grid = pt.grid.Particlepusherinterface(f, 'rho')

