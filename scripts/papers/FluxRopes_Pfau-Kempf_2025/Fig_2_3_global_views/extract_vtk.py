import pytools as pt
import numpy as np
import yt
import matplotlib.pyplot as plt
from evtk.hl import polyLinesToVTK
from multiprocessing import Pool
import pickle
import os

#levels = np.arange(2, 10.1, 2)
levels = (3.0, 5.0, 7.0)
numproc = 4
numthread = 64
chunksize=100

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-t', help="Index of the file to process")
parser.add_argument('-tracing', help="Run tracing (1) or VTK exporting (0)")
args = parser.parse_args()

if args.t == None:
   print("No file index given with the -t option, exiting.")
   exit()
if args.tracing == None:
   print("No value given to -tracing option, exiting.")
   exit()
if args.tracing == "1":
   doTracing = True
else:
   doTracing = False

if doTracing:
   f=pt.vlsvfile.VlsvReader("/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1."+str(args.t).rjust(7, '0')+".vlsv")

   ### Prepare boxes
   [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
   simext=[xmin,xmax,ymin,ymax,zmin,zmax]
   [xsizefs, ysizefs, zsizefs] = f.get_spatial_mesh_size() * 8
   sizesfs = np.array([xsizefs,ysizefs,zsizefs])

   print("sizesfs", sizesfs)

   dx = (xmax-xmin)/xsizefs
   dy = (ymax-ymin)/ysizefs
   dz = (zmax-zmin)/zsizefs

   fullB = f.read_variable_as_fg("vg_b_vol")

   Bx=fullB[:,:,:, 0]
   By=fullB[:,:,:, 1]
   Bz=fullB[:,:,:, 2]
   data=dict(Bx=Bx,By=By,Bz=Bz)

   yt.enable_parallelism()

   yt_dataset = yt.load_uniform_grid(
      data,
      sizesfs,
      bbox=np.array([[xmin - 0.5*dx, xmax - 0.5*dx],
                    [ymin - 0.5*dy, ymax - 0.5*dy],
                    [zmin - 0.5*dz, zmax - 0.5*dz]]),
      periodicity=[True, True, True],
      nprocs=numproc)


   cid=f.read_variable("CellID")
   fluxrope=f.read_variable("vg_fluxrope")
   curvature=f.read_variable("vg_curvature")

   cid_list=cid[np.where(np.logical_and(fluxrope > 0, fluxrope <= levels[-1]))]
   fluxrope_list=fluxrope[np.where(np.logical_and(fluxrope > 0, fluxrope <= levels[-1]))]
   curvature_list=curvature[np.where(np.logical_and(fluxrope > 0, fluxrope <= levels[-1]))]
   array=[]
   for i,val in enumerate(cid_list):
      crd = f.get_cell_coordinates(val)
      c = curvature_list[i]
      c_rad = 1/np.sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2])
      array.append([crd[0],crd[1],crd[2],fluxrope_list[i],c_rad])


   array = np.array(array)

   streamlength = np.max(array[:,3]*array[:,4])

   streamline_seeds = array[:,:3]

   streamlines_pos = yt.visualization.api.Streamlines(yt_dataset, streamline_seeds, "Bx", "By", "Bz", length=streamlength, direction=1, dx=1e6)
   streamlines_neg = yt.visualization.api.Streamlines(yt_dataset, streamline_seeds, "Bx", "By", "Bz", length=streamlength, direction=-1, dx=1e6)

   streamlines_pos.integrate_through_volume();
   streamlines_neg.integrate_through_volume();

   from mpi4py import MPI
   MPI.COMM_WORLD.Barrier()
   if yt.is_root():
      with open("pickledumps/array."+str(args.t)+".pickledump", "wb") as fp:
         pickle.dump(array, fp, protocol=pickle.HIGHEST_PROTOCOL)
      with open("pickledumps/streamlines_pos."+str(args.t)+".pickledump", "wb") as fp:
         pickle.dump(streamlines_pos.streamlines, fp, protocol=pickle.HIGHEST_PROTOCOL)
      with open("pickledumps/streamlines_neg."+str(args.t)+".pickledump", "wb") as fp:
         pickle.dump(streamlines_neg.streamlines, fp, protocol=pickle.HIGHEST_PROTOCOL)
   MPI.COMM_WORLD.Barrier()
else:
   with open("pickledumps/array."+str(args.t)+".pickledump", "rb") as fp:
      array = pickle.load(fp)
   with open("pickledumps/streamlines_pos."+str(args.t)+".pickledump", "rb") as fp:
      streamlines_pos = pickle.load(fp)
   with open("pickledumps/streamlines_neg."+str(args.t)+".pickledump", "rb") as fp:
      streamlines_neg = pickle.load(fp)

   def process_streamlines(indices):
      returnlines=[]
      if np.size(indices) == 1:
         indices=[indices]
      for i in indices:
         for streamlines in (streamlines_pos, streamlines_neg):
            if array[i,3] <= level:
               streamline = streamlines[i]
               streamline = streamline[np.all(streamline != 0.0, axis=1)]
               streamline = np.ndarray.tolist(np.ndarray.flatten(streamline[np.where(np.cumsum(np.linalg.norm(streamline[:-1]-streamline[1:], axis=1)) < 12*array[i,4])]))
               # returnline format [index, coordinates, #points]
               returnlines.append([i, streamline, np.size(streamline)//3])
      if len(returnlines) != 0:
         return returnlines
      else:
         return [[0,0,0]]

   for level in levels:
      if(os.path.exists("streamlines/streamlines_"+str(args.t)+"_"+str(level)+".vtk.vtu")):
         continue

      print("level", level)
            
      print("starting pool")
      ## Parallel processing
      pool = Pool(numthread)
      return_array = pool.imap_unordered(process_streamlines, np.arange(len(streamlines_pos)), chunksize)
#      return_array = pool.map(process_streamlines, np.arange(np.size(streamline_seeds) // 3))
      pool.close()
      pool.join()
      print("stopped pool")

#      return_array = process_streamlines(np.arange(np.size(streamline_seeds) // 3))
      
      storage = []
      curvature = []
      distance = []
      pointsPerLine = []

      for arrayline in return_array:
         if arrayline != [[0,0,0]]:
            for k in np.arange(len(arrayline)):
               for l in arrayline[k][1]:
                  storage.append(l)
               for l in np.arange(arrayline[k][2]):
                  curvature.append(array[arrayline[k][0],4])
                  distance.append(np.sqrt((arrayline[k][1][0]-arrayline[k][1][l*3])**2 + (arrayline[k][1][1]-arrayline[k][1][l*3+1])**2 + (arrayline[k][1][2]-arrayline[k][1][l*3+2])**2))

               pointsPerLine.append(arrayline[k][2])
      
      print("saving VTK")
      storage = np.reshape(np.array(storage), (-1,3))
      # keep np.array() around the storage vectors or it'll complain about stuff not being C and F contiguous.
      polyLinesToVTK("streamlines/streamlines_"+str(args.t)+"_"+str(level)+".vtk", np.array(storage[:,0]), np.array(storage[:,1]), np.array(storage[:,2]), pointsPerLine = np.array(pointsPerLine), pointData = {"curvature_radius" : np.array(curvature), "distance_from_seed" : np.array(distance)})
