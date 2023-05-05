# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

import pytools as pt
import numpy as np
import sys
import matplotlib.pyplot as plt
import scipy
from multiprocessing import shared_memory
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D

# https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])

RE = 6371e3

file_id = 806
path = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/visualizations/lmn/"
fn = "jlsidecar_mva_bulk1.{:07d}.vlsv".format(file_id)
flmn = pt.vlsvfile.VlsvReader(path+fn)
fxo = pt.vlsvfile.VlsvReader(path+"pyXO3/pyXO_bulk1.{:07d}.vlsv".format(file_id))



Ls = flmn.read_variable("LMN_magnetopause/vg_L")
Ms = flmn.read_variable("LMN_magnetopause/vg_M")
Ns = flmn.read_variable("LMN_magnetopause/vg_N")
ds = fxo.read_variable("vg_LMN_neutral_line_distance")
Js = flmn.read_variable("vg_J")
d = flmn.read_variable("LMN_magnetopause/vg_jacobian")
Dimn = fxo.read_variable("vg_MDD_dimensionality")
CellIds = flmn.read_variable("CellID")

picklefn ="/wrk-vakka/group/spacephysics/vlasiator/3D/EGI/visualizations/lmn/806-coords.npz"
try:
   loads = np.load(picklefn)
   B = loads['B']
   dxs = loads['dxs']
   allcoords = loads['allcoords']
   print("Loaded B,dxs")
except:
   B = flmn.read_variable_as_fg("vg_b_vol")
   dxs = flmn.read_variable_as_fg("vg_dxs")
   allcoords = np.array([flmn.get_cell_coordinates(cid) for cid in CellIds])
   np.savez(picklefn, B=B, dxs=dxs,allcoords=allcoords)



outfolder = "/proj/mjalho/analysator/scripts/inspect_xo/"
# dB = pt.calculations.idmesh3d2(CellIds, B, flmn.get_max_refinement_level(), *flmn.get_spatial_mesh_size(), 3)
# ds = pt.calculations.idmesh3d2(CellIds, ds, flmn.get_max_refinement_level(), *flmn.get_spatial_mesh_size(), None)
boxed = np.all(np.abs(allcoords) < 40*6371e3,axis=-1) & (np.sum(allcoords**2,axis=-1) > (5.5*RE)**2)

# these are weird.
sus_cellids = [354515492]

hits = np.array(ds < 0.36) & boxed
scatterpts = allcoords[hits,:]
# print(np.max(dB),np.max(ds))
cid = 355282688
cid = 41099615
for i,cid in enumerate(CellIds): #[355282688]):#,355282689]):
   if(ds[i] > 0.36):
     continue
   #i = np.argwhere(CellIds == cid)
   dx = flmn.get_cell_dx(cid)
   reflevel = flmn.get_amr_level(cid)
   coords = allcoords[i,:]
   if(np.any(np.abs(coords)/6371e3 > 40) or np.linalg.norm(coords) < 5.5*RE):
       continue
   M = (Ms[i,:]/np.linalg.norm(Ms[i,:])).squeeze()
   L = (Ls[i,:]/np.linalg.norm(Ls[i,:])).squeeze()
   N = (Ns[i,:]/np.linalg.norm(Ns[i,:])).squeeze()

   nd = 6
   low = coords - dx*nd
   hi = coords + dx*nd
   cp = np.dot(coords/RE,M)
   # l0 = np.dot(cp - dx*L, L)
   # l1 = np.dot(cp + dx*L, L)
   # n0 = np.dot(cp - dx*N, N)
   # n1 = np.dot(cp + dx*N, N)
   # m0
   # m1
   # print(coords/6371e3, cp, cid, reflevel)
   
   try:
      lowi, upi = flmn.get_bbox_fsgrid_slicemap(low,hi)
      Bsub = flmn.get_bbox_fsgrid_subarray(low,hi,B)
      dsub = dxs[lowi[0]:upi[0]+1, lowi[1]:upi[1]+1, lowi[2]:upi[2]+1, 0]
   except:
      print(cid, "failed at get_bbox_fsgrid_subarray")
      continue   
   l0 = -np.min(dx)*nd
   l1 = +np.min(dx)*nd
   n0 = -np.min(dx)*nd
   n1 = +np.min(dx)*nd
   # print(Ls[i])
   rot = np.stack([L,M,N])
   # print(rot)
   nfg = Bsub[:,:,:,0].shape
   [XmeshXY,YmeshXY] = np.meshgrid(np.linspace(l0/1e3,l1/1e3,num=nfg[0]+1),np.linspace(n0/1e3,n1/1e3,num=nfg[2]+1),indexing='xy')
   XmeshCentres = XmeshXY[:-1,:-1] + 0.5*(XmeshXY[0,1]-XmeshXY[0,0])
   YmeshCentres = YmeshXY[:-1,:-1] + 0.5*(YmeshXY[1,0]-YmeshXY[0,0])

   Bl = np.dot(Bsub,L)
   Bm = np.dot(Bsub,M)
   Bn = np.dot(Bsub,N)
   Bmag = np.linalg.norm(Bsub,axis=-1)

   # Bx = scipy.ndimage.affine_transform(Bsub[:,:,:,0],np.linalg.inv(rot),mode='constant',order=1,offset=(2*nd+1)/2.-np.matmul(np.linalg.inv(rot),np.array([11/2.,11/2.,11/2.])))
   # By = scipy.ndimage.affine_transform(Bsub[:,:,:,1],np.linalg.inv(rot),mode='constant',order=1,offset=(2*nd+1)/2.-np.matmul(np.linalg.inv(rot),np.array([11/2.,11/2.,11/2.])))
   # Bz = scipy.ndimage.affine_transform(Bsub[:,:,:,2],np.linalg.inv(rot),mode='constant',order=1,offset=(2*nd+1)/2.-np.matmul(np.linalg.inv(rot),np.array([11/2.,11/2.,11/2.])))
   # Bxyz = np.stack([Bx,By,Bz],axis=-1)
   #offset = np.array([(2*nd+1)/2.,(2*nd+1)/2.,(2*nd+1)/2.])
   c_in=0.5*np.array(Bl.shape)
   c_out = c_in
   offset=c_in-rot.T@c_out
   Bl = scipy.ndimage.affine_transform(Bl,np.linalg.inv(rot),mode='nearest',order=1,offset=offset)
   Bm = scipy.ndimage.affine_transform(Bm,np.linalg.inv(rot),mode='nearest',order=1,offset=offset)
   Bn = scipy.ndimage.affine_transform(Bn,np.linalg.inv(rot),mode='nearest',order=1,offset=offset)
   # Bmag = scipy.ndimage.affine_transform(Bmag,np.linalg.inv(rot),mode='nearest',order=1,offset=offset)
   Bmag = (Bl**2 + Bn**2)**0.5
   dxsln = scipy.ndimage.affine_transform(dsub,np.linalg.inv(rot),mode='nearest',order=1,offset=offset)
   
   plt.close()

   fig = plt.figure(figsize=(16,12))
   ax_LN_down = fig.add_subplot(4,2,1)
   ax_LN_mid = fig.add_subplot(4,2,3)
   ax_LN_up = fig.add_subplot(4,2,5)
   ax_3d_overview = fig.add_subplot(4,2,(2,4), projection="3d")
   ax_LM = fig.add_subplot(4,2,7)
   ax_3d_detail = fig.add_subplot(4,2,(6,8), projection="3d")
   # try:
   # lic_result = lic.lic(1e9*Bl[:,nfg[1]//2,:],1e9*Bn[:,nfg[1]//2,:], length=20)
   # plt.imshow(lic_result, origin='lower', cmap='gray',extent=(np.min(XmeshXY), np.max(XmeshXY), np.min(YmeshXY), np.max(YmeshXY)))
   # ax1.streamplot(XmeshCentres,YmeshCentres,Bl[:,int(nfg[1]//2+1.1*dx[1]/1e6),:].transpose(),Bn[:,int(nfg[1]//2+1.1*dx[1]/1e6),:].transpose(), linewidth=2,arrowsize=0, color="silver")
   # ax1.streamplot(XmeshCentres,YmeshCentres,Bl[:,int(nfg[1]//2-1.1*dx[1]/1e6),:].transpose(),Bn[:,int(nfg[1]//2-1.1*dx[1]/1e6),:].transpose(), linewidth=2,arrowsize=0, color="gray")
   
   Nslice = (nfg[2]//4,-nfg[2]//4)
   Mslice = int(nfg[1]//2 - dx[1]/1e6)
   ax_LN_down.streamplot(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bl[:,Mslice,Nslice[0]:Nslice[1]].transpose(),Bn[:,Mslice,Nslice[0]:Nslice[1]].transpose())

   pcm = ax_LN_down.pcolormesh(XmeshXY[Nslice[0]:Nslice[1],:],YmeshXY[Nslice[0]:Nslice[1],:], 1e9*Bmag[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="nuuk_r")
   cb = plt.colorbar(pcm,ax=ax_LN_down)
   cb.set_label("$B_{LN}/\mathrm{nT}$")
   ax_LN_down.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],dxsln[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="gray_r", levels=[1.5e6, 3e6, 6e6])
   ax_LN_down.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bl[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="bam", levels=[0.0])
   ax_LN_down.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bn[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="bam_r",levels=[0.0])

   ax_LN_down.set_xlabel("L/km")
   ax_LN_down.set_ylabel("M = {:2.0f} km\nN/km".format(-dx[1]/1e3))
   #ax_LN_down.set_title("M = {:2.0f} km".format(-dx[1]/1e3))
   ax_LN_down.grid()
   cmap = mpl.cm.get_cmap("bam",256)
   custom_lines = [Line2D([0], [0], color="tab:blue", lw=2, label="$B_{LN}$"),
                     Line2D([0], [0], color=cmap(0.0), lw=2, label="$B_L = 0$"),
                     Line2D([0], [0], color=cmap(1.), lw=2, label="$B_N = 0$"),
                     Line2D([0], [0], color="gray", lw=2,label="$\Delta{}X_{vg}$")]
   ax_LN_down.legend(handles=custom_lines,loc="upper right")


   Mslice = int(nfg[1]//2)
   ax_LN_mid.streamplot(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bl[:,Mslice,Nslice[0]:Nslice[1]].transpose(),Bn[:,Mslice,Nslice[0]:Nslice[1]].transpose())

   pcm = ax_LN_mid.pcolormesh(XmeshXY[Nslice[0]:Nslice[1],:],YmeshXY[Nslice[0]:Nslice[1],:], 1e9*Bmag[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="nuuk_r")
   cb = plt.colorbar(pcm,ax=ax_LN_mid)
   cb.set_label("$B_{LN}/\mathrm{nT}$")
   ax_LN_mid.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],dxsln[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="gray_r", levels=[1.5e6, 3e6, 6e6])
   ax_LN_mid.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bl[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="bam", levels=[0.0])
   ax_LN_mid.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bn[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="bam_r",levels=[0.0])

   ax_LN_mid.set_xlabel("L/km")
   ax_LN_mid.set_ylabel("M = 0 km\nN/km")
   #ax_LN_mid.set_title("M = 0 km")
   ax_LN_mid.legend(handles=custom_lines,loc="upper right")
   ax_LN_mid.grid()

   Mslice = int(nfg[1]//2 + dx[1]/1e6)
   ax_LN_up.streamplot(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bl[:,Mslice,Nslice[0]:Nslice[1]].transpose(),Bn[:,Mslice,Nslice[0]:Nslice[1]].transpose())

   pcm = ax_LN_up.pcolormesh(XmeshXY[Nslice[0]:Nslice[1],:],YmeshXY[Nslice[0]:Nslice[1],:], 1e9*Bmag[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="nuuk_r")
   cb = plt.colorbar(pcm,ax=ax_LN_up)
   cb.set_label("$B_{LN}/\mathrm{nT}$")
   ax_LN_up.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],dxsln[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="gray_r", levels=[1.5e6, 3e6, 6e6])
   ax_LN_up.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bl[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="bam", levels=[0.0])
   ax_LN_up.contour(XmeshCentres[Nslice[0]:Nslice[1],:],YmeshCentres[Nslice[0]:Nslice[1],:],Bn[:,Mslice,Nslice[0]:Nslice[1]].transpose(),cmap="bam_r",levels=[0.0])

   ax_LN_up.set_xlabel("L/km")
   ax_LN_up.set_ylabel("M = {:2.0f} km\nN/km".format(dx[1]/1e3))
   #ax_LN_up.set_title("M = {:2.0f} km".format(dx[1]/1e3))
   ax_LN_up.grid()
   custom_lines = [Line2D([0], [0], color="tab:blue", lw=2, label="$B_{LN}$"),
                     Line2D([0], [0], color=cmap(0.0), lw=2, label="$B_L = 0$"),
                     Line2D([0], [0], color=cmap(1.), lw=2, label="$B_N = 0$"),
                     Line2D([0], [0], color="gray", lw=2,label="$\Delta{}X_{vg}$")]
   ax_LN_up.legend(handles=custom_lines,loc="upper right")

   # except:
   #    print(cid, "failed at plotting")
   #    continue

   ax_LM.streamplot(XmeshCentres,YmeshCentres,Bl[:,:,nfg[2]//2].transpose(),Bm[:,:,nfg[2]//2].transpose())
   pcm = ax_LM.pcolormesh(XmeshXY,YmeshXY, 1e9*Bm[:,:,nfg[2]//2].transpose(),norm=mpl.colors.CenteredNorm(),cmap="vik")
   cb = plt.colorbar(pcm,ax=ax_LM)
   cb.set_label("$B_M/\mathrm{nT}$")
   ax_LM.contour(XmeshCentres,YmeshCentres,dxsln[:,:,nfg[2]//2].transpose(),cmap="gray_r", levels=[1.5e6, 3e6, 6e6])
   ax_LM.contour(XmeshCentres,YmeshCentres,Bl[:,:,nfg[2]//2].transpose(),cmap="bam", levels=[0.0])
   ax_LM.contour(XmeshCentres,YmeshCentres,Bn[:,:,nfg[2]//2].transpose(),cmap="bam_r",levels=[0.0])
   ax_LM.set_xlabel("L/km")
   ax_LM.set_ylabel("N = 0 km\nM/km")
   #ax_LM.set_title("N = 0 km")
   ax_LM.grid()
   ax_LM.legend(handles=custom_lines,loc="upper right")


   u = np.linspace(0, 2 * np.pi, 100)
   v = np.linspace(0, np.pi, 100)
   x = 1 * np.outer(np.cos(u), np.sin(v))
   y = 1 * np.outer(np.sin(u), np.sin(v))
   z = 1 * np.outer(np.ones(np.size(u)), np.cos(v))
   ax_3d_overview.plot_surface(x, y, z)
   sc = ax_3d_overview.scatter(allcoords[hits,0]/RE,allcoords[hits,1]/RE,allcoords[hits,2]/RE, c=d[hits,6], s=1, alpha=0.3, marker=',', cmap="roma_r",vmin=-1e-6,vmax=1e-6,label="X/O")
   cb = plt.colorbar(sc,ax=ax_3d_overview, pad=0.2)
   cb.set_label("$\partial{}_LB_N / \mathrm{Tm}^-1$ (positive is X)")
   #ax2.set_aspect('equal')
   ax_3d_overview.quiver(coords[0]/RE,coords[1]/RE,coords[2]/RE, L.squeeze()[0], L.squeeze()[1], L.squeeze()[2], color="r", length=20, label="L")
   ax_3d_overview.quiver(coords[0]/RE,coords[1]/RE,coords[2]/RE, M.squeeze()[0], M.squeeze()[1], M.squeeze()[2], color="g", length=20, label="M")
   ax_3d_overview.quiver(coords[0]/RE,coords[1]/RE,coords[2]/RE, N.squeeze()[0], N.squeeze()[1], N.squeeze()[2], color="b", length=20, label="N")
   ax_3d_overview.set_xlabel("$X/\mathrm{R_E}$")
   ax_3d_overview.set_ylabel("$Y/\mathrm{R_E}$")
   ax_3d_overview.set_zlabel("$Z/\mathrm{R_E}$")
   ax_3d_overview.legend(loc="upper right")
   set_axes_equal(ax_3d_overview)

   hits_small_x = hits & np.all(np.abs(allcoords-coords) < dx*nd,axis=-1) & (np.sum(allcoords**2,axis=-1) > (5.5*RE)**2) & (d[:,6] > 0)
   hits_small_o = hits & np.all(np.abs(allcoords-coords) < dx*nd,axis=-1) & (np.sum(allcoords**2,axis=-1) > (5.5*RE)**2) & (d[:,6] < 0)

   ncoords = (rot@((allcoords-coords).T)).T

   sc = ax_3d_detail.scatter(ncoords[hits_small_x,0]/RE,ncoords[hits_small_x,1]/RE,ncoords[hits_small_x,2]/RE, c=d[hits_small_x,6], s=40, marker='x', cmap="roma_r",norm=mpl.colors.CenteredNorm(),label="X")
   sc = ax_3d_detail.scatter(ncoords[hits_small_o,0]/RE,ncoords[hits_small_o,1]/RE,ncoords[hits_small_o,2]/RE, c=d[hits_small_o,6], s=10, marker='o', cmap="roma_r",norm=mpl.colors.CenteredNorm(),label="O")
   cb = plt.colorbar(sc,ax=ax_3d_detail, pad=0.2)
   cb.set_label("$\partial{}_LB_N / \mathrm{Tm}^-1$ (positive is X)")
   #ax2.set_aspect('equal')
   # ax_3d_detail.quiver(coords[0]/RE,coords[1]/RE,coords[2]/RE, L.squeeze()[0], L.squeeze()[1], L.squeeze()[2], color="r", length=3, label="L")
   # ax_3d_detail.quiver(coords[0]/RE,coords[1]/RE,coords[2]/RE, M.squeeze()[0], M.squeeze()[1], M.squeeze()[2], color="g", length=3, label="M")
   # ax_3d_detail.quiver(coords[0]/RE,coords[1]/RE,coords[2]/RE, N.squeeze()[0], N.squeeze()[1], N.squeeze()[2], color="b", length=3, label="N")
   Ll = rot@L
   Mm = rot@M
   Nn = rot@N
   ax_3d_detail.quiver(0,0,0, Ll.squeeze()[0], Ll.squeeze()[1], Ll.squeeze()[2], color="r", length=3, label="L")
   ax_3d_detail.quiver(0,0,0, Mm.squeeze()[0], Mm.squeeze()[1], Mm.squeeze()[2], color="g", length=3, label="M")
   ax_3d_detail.quiver(0,0,0, Nn.squeeze()[0], Nn.squeeze()[1], Nn.squeeze()[2], color="b", length=3, label="N")
   ax_3d_detail.set_xlabel("$L/\mathrm{R_E}$")
   ax_3d_detail.set_ylabel("$M/\mathrm{R_E}$")
   ax_3d_detail.set_zlabel("$N/\mathrm{R_E}$")
   ax_3d_detail.legend(loc="upper right")
   set_axes_equal(ax_3d_detail)

   xostr = ""
   if (d[i,6].squeeze()> 0):
      xostr = " (X)"
      if (d[i,6].squeeze()*dx[0] > 1e-6*1e6):
         folder = "X_1e-0"
      elif (d[i,6].squeeze()*dx[0] > 1e-7*1e6):
         folder = "X_1e-1"
      elif (d[i,6].squeeze()*dx[0] > 1e-8*1e6):
         folder = "X_1e-2"
      elif (d[i,6].squeeze()*dx[0] > 1e-9*1e6):
         folder = "X_1e-3"
      else:
         folder = "X_low"
   else:
      xostr = " (O)"
      if (d[i,6].squeeze()*dx[0] < -1e-6*1e6):
         folder = "O_1e-0"
      elif (d[i,6].squeeze()*dx[0] < -1e-7*1e6):
         folder = "O_1e-1"
      elif (d[i,6].squeeze()*dx[0] < -1e-8*1e6):
         folder = "O_1e-2"
      elif (d[i,6].squeeze()*dx[0] < -1e-9*1e6):
         folder = "O_1e-3"
      else:
         folder = "O_low"

   #plt.title(["r = {0:3f}".format(*tuple(coords/RE)), "\n (LNM) = (" + str(L) +", "+ str(M) + "," +str(N)])
   plt.suptitle("CellID {:d}".format(cid) + ", J = ({:3e},{:3e},{:3e}), {:3e} A".format(*tuple(Js[i,:]),np.linalg.norm(Js[i,:])) +
               "\nr/RE = [{0:.1f}, {1:.1f}, {2:.1f}]; dBNdL = {3:.3e}".format(*tuple(coords/RE),d[i,6].squeeze()) + xostr + ", ds = {:.2f}".format(ds[i]) + ", Di = ({0:.2f},{1:.2f},{2:.2f})".format(*tuple(Dimn[i,:])) +
               "\n LMN: ({0:.2f},{1:.2f},{2:.2f}),({3:.2f},{4:.2f},{5:.2f}),({6:.2f},{7:.2f},{8:.2f})".format(*tuple(L),*tuple(M),*tuple(N)))
   
   plt.savefig("806/"+folder+"/fig{:012d}.png".format(cid))
   # print(rot @ np.array([[1,0,0],[0,1,0],[0,0,1]]))

# print(lowi,upi)

# pt.plot.plot_colormap3dslice(path+fn,
#                               outputdir=outfolder,
#                               var='vg_J',
#                               streamlines='vg_b_vol',
#                               normal = M,#N,
#                               up = N,
#                               useimshow=True,
#                               cutpointre = cp,
#                               #boxm=[x0, x1, y0, y1 ] )
   #   fig1 = ax1.imshow(datamap,
   #                     cmap=colormap,
   #                     norm=norm,
   #                     interpolation=imshowinterp,
   #                     origin='lower',
   #                     extent=(np.min(XmeshPass), np.max(XmeshPass), np.min(YmeshPass), np.max(YmeshPass))
   #                    )
                              #streamplot
   
