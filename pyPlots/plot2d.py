# 2D plotting

def vlsv_plot2d(vlsvReader,varName="rho",unit=1,xFigSize=10,yFigSize=10):
 ''' 2-dimensional pseudocolor plot
 :param vlsvReader: VLSV reader interface
 :param varName:    Variable name if VLSV reader
 :param unit:       Length unit of axes
 :param xFigSize:   Figure size in horizontal direction
 :param yFigSize:   Figure size in vertical direction
 :returns: nothing
 .. code-block:: python
 # Example usage:
   import pytools as pt
   fvlsv = pt.vlsvfile.VlsvReader("/path/to/run/bulk.0003000.vlsv")
   pt.plot.vlsv_plot2d(fvlsv,varName="rho",unit=6371e3)
  '''
 import numpy as np
 import matplotlib.pyplot as plt
 import operator as oper
 # read file header
 xmin = vlsvReader.read_parameter("xmin")
 xmax = vlsvReader.read_parameter("xmax")
 ymin = vlsvReader.read_parameter("ymin")
 ymax = vlsvReader.read_parameter("ymax")
 zmin = vlsvReader.read_parameter("zmin")
 zmax = vlsvReader.read_parameter("zmax")
 xsize = xmax - xmin
 ysize = ymax - ymin
 zsize = zmax - zmin
 nx = vlsvReader.read_parameter("xcells_ini")
 ny = vlsvReader.read_parameter("ycells_ini")
 nz = vlsvReader.read_parameter("zcells_ini")
 dx = xsize/nx
 dy = ysize/ny
 dz = zsize/nz
 # check if the file has xy or xz plane
 if(ny > nz and nz == 1):
  plane = "xy"
  meshX,meshY = np.meshgrid(np.linspace(xmin/unit,xmax/unit,nx),np.linspace(ymin/unit,ymax/unit,ny))
 elif(ny < nz and ny == 1):
  plane = "xz"
  meshX,meshY = np.meshgrid(np.linspace(xmin/unit,xmax/unit,nx),np.linspace(zmin/unit,zmax/unit,nz))
 else:
  print "Error: cannot determine simulation plane"
  return
 # cell id - index  dict
 locs = vlsvReader.get_cellid_locations()
 cellids = locs.keys()
 # open variable
 if varName in vlsvReader.get_all_variables():
  var = vlsvReader.read_variable(varName)
 else:
  print "Variable " + varName + " not found"
  return
 # sort variable array according to cell ids
 locs_sorted = sorted(locs.iteritems(), key=oper.itemgetter(0))
 var_sorted = []
 for i in locs_sorted:
  var_sorted.append(var[i[1]])
 # Shape variable array as 2-dimensional matrix
 if(plane == "xz"):
  var = np.asarray(var_sorted).reshape(nx,nz) 
 if(plane == "xy"):
  var = np.asarray(var_sorted).reshape(nx,ny)
 # do plotting
 plt.figure(figsize=(xFigSize,yFigSize))
 plt.pcolormesh(meshX,meshY,var)
 plt.axis('equal')
 plt.xlabel("x [" + str(unit/1e3) + " km]")
 plt.xlim([xmin/unit,xmax/unit])
 if(plane == "xz"):
  plt.ylim([zmin/unit,zmax/unit])
  plt.ylabel("z [" + str(unit/1e3) + " km]")
 if(plane == "xy"):
  plt.ylim([ymin/unit,ymax/unit])
  plt.xlabel("y [" + str(unit/1e3) + " km]")
 plt.show()

