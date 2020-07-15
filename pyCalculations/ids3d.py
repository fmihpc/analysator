import numpy as np


# finds the cell ids which are needed to plot a 2d cut through out of a 3d mesh # WARNING! works only if cellids is sorted
def ids3d(cellids, depth, reflevel,
          xsize, ysize, zsize,
          xmin=None, xmax=None,
          ymin=None, ymax=None,
          zmin=None, zmax=None):

  # is the cut in x, y or z direction
  if xmin != None and xmax != None:
    imin = xmin; imax = xmax; size = xsize; sett = 0
  elif ymin != None and ymax != None:
    imin = ymin; imax = ymax; size = ysize; sett = 1
  elif zmin != None and zmax != None:
    imin = zmin; imax = zmax; size = zsize; sett = 2

  # find the cut through depths (index) for each refinement level
  pro = depth/(imax-imin)
  depths = []
  for i in range(reflevel+1):
    depth = int(pro*size*2**i) + 1 # goes from 1 to Nmax
    if depth > int(size*2**i):
      print("depth error ",depth,i)
      depth -= 1
    depths.append(depth)

  #############
  # find the ids
  length = 0
  cellids = np.array(cellids)
  idlist = np.array([]); indexlist = np.array([])
  cells = int(xsize*ysize*zsize); cellsum = cells; cellsumII = 0 # cellsum = (the number of cells up to refinement 
  # level i); cellsumII = (the number of cells up to refinement level i-1) 
  for i in range(reflevel+1):
    # cell ids at the refinement level i
    ids = cellids[(cellsumII < cellids) & (cellids <= cellsum)]

    # compute every cell ids x, y and z coordinates
    z = (ids - cellsumII - 1)//(xsize*ysize*4**i) +1 # goes from 1 to NZ
    #z = np.where(np.isclose(z, z.astype(int)), z, z + 1)
    z = z.astype(int)

    # cellsumB = (cellsumII + the number of ids up to the coordinate z in the refinement level i)
    cellsumB = ((z-1)*xsize*ysize*4**i).astype(int) + cellsumII

    y = (ids - cellsumB - 1)//(xsize*2**i) +1 # goes from 1 to NY
    #y = np.where(y - y.astype(int) == 0, y, y + 1)
    y = y.astype(int)

    cellsumC = ((y-1)*xsize*2**i)
    x = (ids - cellsumB - cellsumC).astype(int) # goes from 1 to NX

    # checks is the cut throughs normal in x, y or z direction
    if sett == 0:
      xyz = x
    elif sett == 1:
      xyz = y
    elif sett == 2:
      xyz = z
      
    # finds the needed elements to create the asked cut throuhg and puts the results in the indexlist and the idlist
    elements = xyz == depths[i]
    indexlist = np.append(indexlist, np.arange(length, length+len(xyz))[elements])
    idlist = np.append(idlist, ids[elements])

    # update these values to match refinement level i+1
    length += len(xyz)
    cellsumII = cellsum
    cellsum += cells*8**(i+1)
  # returns a list of ids (idlist) and a equivalent list of indices (indexlist)
  return idlist.astype(int), indexlist.astype(int)



# creates 2d grid for ploting
def idmesh3d(idlist, data, reflevel, xsize, ysize, zsize, xyz, datadimension):

  # is the cut in x, y or z direction
  if xyz == 0:
    dsize = int(ysize*zsize*4**reflevel)
    asize = ysize; bsize = zsize
  elif xyz == 1:
    dsize = int(xsize*zsize*4**reflevel)
    asize = xsize; bsize = zsize
  elif xyz == 2:
    dsize = int(xsize*ysize*4**reflevel)
    asize = xsize; bsize = ysize

  # datadimension is None for scalar,
  # N for vector, and (N,M) for tensor data
  if datadimension is None:
    dpoints = np.zeros(dsize)
  elif np.ndim(datadimension) == 0:
    dpoints = np.zeros((dsize,datadimension))
  elif np.ndim(datadimension) == 1:
    dpoints = np.zeros((dsize,datadimension[0],datadimension[1]))
  else:
    print("Error finding data dimension in idmesh3d")
    return -1

  ######################
  # create the plot grid
  cells = int(xsize*ysize*zsize)
  cellsum = cells # the number of cells up to refinement level i
  cellsumII = 0 # the number of cells up to refinement level i-1 
  for i in range(reflevel+1):
    # the cell ids and the data at the refinement level i
    ids = idlist[(cellsumII < idlist) & (idlist <= cellsum)]
    if datadimension is None:
      dat = data[(cellsumII < idlist) & (idlist <= cellsum)]
    elif np.ndim(datadimension) == 0:
      dat = data[(cellsumII < idlist) & (idlist <= cellsum),:]
    elif np.ndim(datadimension) == 1:
      dat = data[(cellsumII < idlist) & (idlist <= cellsum),:,:]

    # compute every cell ids x, y and z coordinates
    z = (ids - cellsumII - 1)//(xsize*ysize*4**i) +1 # goes from 1 to NZ
    #z = np.where(np.isclose(z, z.astype(int)), z, z + 1)
    z = z.astype(int)

    # cellsumB = (cellsumII + the number of ids up to the coordinate z in the refinement level i)
    cellsumB = ((z-1)*xsize*ysize*4**i).astype(int) + cellsumII

    y = (ids - cellsumB - 1)//(xsize*2**i) +1 # goes from 1 to NY
    #y = np.where(y - y.astype(int) == 0, y, y + 1)
    y = y.astype(int)

    cellsumC = ((y-1)*xsize*2**i)
    x = (ids - cellsumB - cellsumC).astype(int) # goes from 1 to NX
    
    # gets the correct coordinate values and the widths for the plot
    if xyz == 0:
      a = y; b = z
      l = ysize*2**i
    elif xyz == 1:
      a = x; b = z
      l = xsize*2**i
    elif xyz == 2:
      a = x; b = y
      l = xsize*2**i

    # computes what are the coordinate values in numpy 1d array
    idc = ((a-1 + l*(b-1)*2**(reflevel-i))*2**(reflevel-i)).astype(int)

    # inserts the data values into dpoints
    for j in range(2**(reflevel-i)):
      for k in range(2**(reflevel-i)):
        if datadimension is None:
          dpoints[(idc + j*asize*2**reflevel + k).astype(int)] = dat
        elif np.ndim(datadimension) == 0:
          dpoints[(idc + j*asize*2**reflevel + k).astype(int),:] = dat
        elif np.ndim(datadimension) == 1:
          dpoints[(idc + j*asize*2**reflevel + k).astype(int),:,:] = dat

    # update these values to match refinement level i+1
    cellsumII = cellsum
    cellsum += cells*8**(i+1)

  # reshape dpoints to 2d grid
  if datadimension is None:
    dpoints = dpoints.reshape(int(bsize*2**reflevel), int(asize*2**reflevel))
  elif np.ndim(datadimension) == 0:
    dpoints = dpoints.reshape(int(bsize*2**reflevel), int(asize*2**reflevel),datadimension)
  elif np.ndim(datadimension) == 1:
    dpoints = dpoints.reshape(int(bsize*2**reflevel), int(asize*2**reflevel),datadimension[0],datadimension[1])
      
  # returns a data grid (dpoints) 
  return dpoints


# find the highest refinement level 
def refinement_level(xsize, ysize, zsize, bigid):
  cells = int(xsize*ysize*zsize)
  reflevel = 0
  cellsum = cells
  while bigid > cellsum:
    reflevel += 1
    cellsum += int(cells*8**reflevel)

  # return the highest refinement level
  return reflevel
