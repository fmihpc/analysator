





import analysator as pt
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import uniform_filter1d
import os
from analysator.calculations.lineout import lineout

r_e = 6.371e6 #Maybe we should add an import that holds all the constants?

def jplots(
    var,
    fnr1,
    fnr2,
    bulkpath,
    bulkprefix,
    point1,
    point2,
    outputname,
    outputdir=None,
    filt=-1,
    op="pass",
    cmap="viridis",
    re=False,
    npoints=100,
    interpolation_order=1,
    draw=False,
    nooverwrite=False
):
    '''
    Function for plotting a cut-through timeseries, aka keogram, aka time-elongation map, etc. for a given variable, a given set of coordinates, and a range of times.
    
    :kwarg var: Variable to plot (e.g. "proton/vg_rho", "vg_b_vol")
    :kwarg fnr1: First file number to plot
    :kwarg fnr2: Last file number to plot
    :kwarg bulkpath: Path to directory containing bulk files
    :kwarg bulkprefix: Starting string of bulk file name (e.g. bulk, bulk1, bulk5)
    :kwarg point1: First point on the line to plot
    :kwarg point2: Last point on the line to plot
    :kwarg outputname: Name of output file
    :kwarg outputdir: Path to output directory
    :kwarg filt: Filter out temporally slowly changing signal? (<=0: no filtering, >0: filter with specified window size), default=-1
    :kwarg op: Variable operator, default="pass"
    :kwarg cmap: Colormap, default="viridis"
    :kwarg re: Input points given in Earth radii? default=False
    :kwarg npoints: Number of points in line, default=100
    :kwarg interpolation_order: Order of interpolation (0 or 1), default=1
    :kward nooverwrite: Whether to overwrite if output file already exists, default=False
    :kwarg draw: Whether draw on screen or output to file, default=False

    :returns:           Outputs an image to a file or to the screen.

    .. code-block:: python

        # Example usage:
        jplots(var="vg_b_vol",fnr1=600,fnr2=650,bulkpath='/wrk-vakka/group/spacephysics/vlasiator/3D/FIF/bulk1',bulkprefix='bulk1',
                point1=[10,-20,0],point2=[20,-20,0],outputdir='./',outputname='test.png',op='x',re=True,npoints=50)
    '''


    fnr_arr = np.arange(fnr1, fnr2 + 0.1, 1, dtype=int)
    t_arr = np.zeros_like(fnr_arr).astype(float)

    fobj = pt.vlsvfile.VlsvReader(
        os.path.join(bulkpath , bulkprefix + ".{}.vlsv".format(str(fnr1).zfill(7)))
    )
    if re:
        point1 = [point1[0] * r_e, point1[1] * r_e, point1[2] * r_e]
        point2 = [point2[0] * r_e, point2[1] * r_e, point2[2] * r_e]

    lineout0 = lineout(
        fobj,
        point1,
        point2,
        variable=var,
        operator=op,
        points=npoints,
        interpolation_order=interpolation_order,
    )
    distances, _, _ = lineout0

    fobj.optimize_clear_fileindex_for_cellid()

    distances = np.array(distances) / r_e

    data_arr = np.zeros((fnr_arr.size, distances.size), dtype=float)

    for idx in range(fnr_arr.size):
        fnr = fnr_arr[idx]
        vlsvobj = pt.vlsvfile.VlsvReader(
            os.path.join(bulkpath , bulkprefix + ".{}.vlsv".format(str(fnr).zfill(7)))
        )
        t_arr[idx] = vlsvobj.read_parameter("time")

        linecut = lineout(
            vlsvobj,
            point1,
            point2,
            variable=var,
            operator=op,
            points=npoints,
            interpolation_order=interpolation_order,
        )
        data_arr[idx, :] = linecut[2]
        vlsvobj.optimize_clear_fileindex_for_cellid()
    if filt > 0:
        data_arr = data_arr - uniform_filter1d(data_arr, size=filt, axis=0)

    XmeshXY, YmeshXY = np.meshgrid(distances, t_arr)

    fig, ax = plt.subplots(1, 1, figsize=(8, 12), constrained_layout=True)

    im = ax.pcolormesh(
        XmeshXY,
        YmeshXY,
        data_arr,
        shading="gouraud",
        cmap=cmap,
        rasterized=True,
    )

    ax.set_xlim(distances[0], distances[-1])
    ax.set_ylim(t_arr[0], t_arr[-1])
    ax.set_xlabel("Distance along cut [RE]", labelpad=10, fontsize=16)
    ax.set_ylabel("Time [s]", labelpad=10, fontsize=16)
    ax.set_title(var, pad=10, fontsize=16)

    cb = fig.colorbar(im, ax=ax)

    if not draw:
        outputpath=pt.plot.output_path(outputname,None,outputdir,nooverwrite)
        fig.savefig(outputpath, dpi=300)
        plt.close(fig)
    else:
        return (fig,ax,XmeshXY,YmeshXY,data_arr)


