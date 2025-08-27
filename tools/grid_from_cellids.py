#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

base_grid = np.array((4,4,4))
base_cell_count = base_grid[0] * base_grid[1] * base_grid[2]

explode_gap = 0.1

alpha=0.3
color_cycler = ((1,0,0,alpha), (0,1,0,alpha), (0,0,1,alpha), (1,1,0,alpha), (1,0,1,alpha), (0,1,1,alpha))

def get_amr_level(cellid):
   '''Returns the AMR level of a given cell defined by its cellid

   :param cellid:        The cell's cellid
   :returns:             The cell's refinement level in the AMR
   '''
   stack = True
   if not hasattr(cellid,"__len__"):
      cellid = np.atleast_1d(cellid)
      stack = False

   AMR_count = np.zeros(np.array(cellid).shape, dtype=np.int64)
   cellids = cellid.astype(np.int64)
   iters = 0
   while np.any(cellids > 0):
      mask = cellids > 0
      sub = 2**(3*AMR_count)*base_cell_count
      np.subtract(cellids, sub.astype(np.int64), out = cellids, where = mask)
      np.add(AMR_count, 1, out = AMR_count, where = mask)
      iters = iters+1
      if(iters > 10):
         print("Can't have 10 levels of refinement, can we?")
         break

   if stack:
      return AMR_count - 1
   else:
      return (AMR_count - 1)[0]


def get_cell_indices(cellids, reflevels=None):
   ''' Returns a given cell's indices as a numpy array

   :param cellid:            The cell's ID, numpy array
   :param reflevel:          The cell's refinement level in the AMR
   :returns: a numpy array with the coordinates

   .. note:: The cell ids go from 1 .. max not from 0
   '''

   stack = True
   if not hasattr(cellids,"__len__"):
      cellids = np.atleast_1d(cellids)
      stack = False

   if reflevels is None:
      reflevels = get_amr_level(cellids)
   else:
      reflevels = np.atleast_1d(reflevels)

   mask = reflevels >= 0

   max_refinement_level = np.max(reflevels)

   # Calculating the index of the first cell at this reflevel
   index_at_reflevel = np.zeros(max_refinement_level+1, dtype=np.int64)
   isum = 0
   for i in range(0,max_refinement_level):
      isum = isum + 2**(3*i) * base_cell_count
      index_at_reflevel[i+1] = isum

   # Get cell indices:
   cellids = np.array(cellids - 1 - index_at_reflevel[reflevels], dtype=np.int64)
   cellindices = np.full((len(cellids),3), -1)
   cellindices[mask,0] = (cellids[mask])%(np.power(2,reflevels[mask])*base_grid[0])
   cellindices[mask,1] = ((cellids[mask])//(np.power(2,reflevels[mask])*base_grid[0]))%(np.power(2,reflevels[mask])*base_grid[1])
   cellindices[mask,2] = (cellids[mask])//(np.power(4,reflevels[mask])*base_grid[0]*base_grid[1])

   # Return the indices:
   if stack:
      return np.array(cellindices)
   else:
      return np.array(cellindices)[0]


def explode(data):
   size = np.array(data.shape)*2
   data_e = np.zeros(size - 1, dtype=data.dtype)
   data_e[::2, ::2, ::2] = data
   return data_e



def main():
   try:
      cellids = np.loadtxt("cellids.txt")
   except:
      print("No cellids.txt found, using default cellids.")
      cellids = np.asarray([2,3,436,437,438,6345,6346])


   ref_levels = get_amr_level(cellids)
   max_ref_level = np.max(ref_levels)
   all_indices = get_cell_indices(cellids, ref_levels)


   fig = plt.figure()
   ax = fig.add_subplot(projection='3d')

   legend_handles=[]

   for ref_level in range(max_ref_level+1):
      x, y, z = np.indices(base_grid*2**ref_level)
      #xx, yy, zz = np.indices(base_grid*2**ref_level + 1) / 2**ref_level
      this_level_indices = all_indices[np.argwhere(ref_levels == ref_level)]
      this_level_voxels = np.full_like(x, False)
      for indices in this_level_indices:
         this_level_voxels |= (x == indices[0][0]) & (y == indices[0][1]) & (z == indices[0][2])

      # upscale the above voxel image, leaving gaps
      this_level_voxels_exploded = explode(this_level_voxels)

      # Shrink the gaps
      xx, yy, zz = ((np.indices(np.array(this_level_voxels_exploded.shape) + 1).astype(float) ) // 2 )
      xx[0::2, :, :] += explode_gap / 2**ref_level
      yy[:, 0::2, :] += explode_gap / 2**ref_level
      zz[:, :, 0::2] += explode_gap / 2**ref_level
      xx[1::2, :, :] += 1.0 - explode_gap / 2**ref_level
      yy[:, 1::2, :] += 1.0 - explode_gap / 2**ref_level
      zz[:, :, 1::2] += 1.0 - explode_gap / 2**ref_level

      xx /= 2**ref_level
      yy /= 2**ref_level
      zz /= 2**ref_level

      vx = ax.voxels(xx, yy, zz, this_level_voxels_exploded, edgecolor='k', facecolor=color_cycler[ref_level % len(color_cycler)])
      the_patch = mpatches.Patch(facecolor=color_cycler[ref_level % len(color_cycler)], edgecolor='k', label="refinement level "+str(ref_level))
      legend_handles.append(the_patch)

   ax.legend(handles=legend_handles)

   ax.set_box_aspect((np.ptp(xx), np.ptp(yy), np.ptp(zz)))

   plt.show()


if __name__ == "__main__":
   main()
