"""
An example of using find_find_magnetopause_3d_sw_streamline() and plotting the output.
"""

import numpy as np
import matplotlib.pyplot as plt
import analysator as pt

## get bulk data
run =  'EGI' # for plot output file name
file_name="/wrk-vakka/group/spacephysics/vlasiator/3D/"+run+"/bulk/bulk5.0000070.vlsv" 

## magnetopause
x_point_n = 50 # number needed for 2d slice plots, default is 50

vertices, faces = pt.calculations.find_magnetopause_sw_streamline_3d(
    file_name, 
    streamline_seeds=None,
    seeds_n=25,
    seeds_x0=20*6371000, 
    seeds_range=[-5*6371000, 5*6371000],
    dl=2e6, 
    iterations=200, 
    end_x=-15*6371000,
    x_point_n=x_point_n,
    sector_n=36)


outdir="" # output plot save directory

### 3D surface plot ###

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_trisurf(vertices[:, 0], vertices[:,1], faces, vertices[:,2], linewidth=0.2)
plt.savefig(outdir+run+'_magnetopause_3D_surface.png')
plt.close()


### 2D slice plots ###

vertices = vertices/6371000 # to RE
magnetopause = np.array(np.split(vertices, x_point_n)) # magnetopause points grouped by x-axis coordinate

## take separate arrays for different 2d slice plots

quarter_slice = int(np.shape(magnetopause)[1]/4)
# xy plane: z=0
xy_slice = np.concatenate((magnetopause[:,0][::-1], magnetopause[:,2*quarter_slice]))
# xz plane: y=0
xz_slice = np.concatenate((magnetopause[:,quarter_slice][::-1], magnetopause[:,3*quarter_slice]))


# analysator 3dcolormapslice y=0

def external_plot(ax,XmeshXY=None, YmeshXY=None, pass_maps=None, requestvariables=False):
    if requestvariables==True:
        return ['vg_v']
    ax.plot(xz_slice[:,0], xz_slice[:,2],  color='limegreen', linewidth=1.5)

pt.plot.plot_colormap3dslice(
filename=file_name,
run=run,
outputdir=outdir,
boxre = [-21,21,-21,21],
external = external_plot,
colormap='inferno',
normal='y'
)

# analysator 3dcolormapslice z=0

def external_plot(ax,XmeshXY=None, YmeshXY=None, pass_maps=None, requestvariables=False):
    if requestvariables==True:
        return ['vg_v']
    ax.plot(xy_slice[:,0], xy_slice[:,1], color='limegreen', linewidth=1.5)

pt.plot.plot_colormap3dslice(
filename=file_name,
outputdir=outdir,
run=run,
boxre = [-21,21,-21,21],
external = external_plot,
colormap='inferno',
normal='z')
