import sys, os
import numpy as np
from scipy.ndimage import uniform_filter1d
import pytools as pt
import matplotlib.pyplot as plt
import argparse

r_e = 6.371e6


def jplots(
    var,
    fnr0,
    fnr1,
    start_coords,
    end_coords,
    dr,
    bulkpath,
    bulkprefix,
    outputname,
    outputdir,
    intpol=False,
    filt=-1,
    op="pass",
    cmap="viridis",
    pointfile=None,
):

    # dr *= 1000
    # dr /= r_e

    fnr_arr = np.arange(fnr0, fnr1 + 0.1, 1, dtype=int)
    t_arr = np.zeros_like(fnr_arr).astype(float)

    if bulkpath[-1] != "/":
        bulkpath += "/"

    if pointfile is None:
        x0, y0, z0 = start_coords * r_e / 1000
        x1, y1, z1 = end_coords * r_e / 1000
        npoints = (
            int(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2) / dr) + 1
        )

        xlist = np.linspace(x0, x1, npoints, dtype=float)
        ylist = np.linspace(y0, y1, npoints, dtype=float)
        zlist = np.linspace(z0, z1, npoints, dtype=float)

        coords = r_e * np.array([xlist, ylist, zlist]).T
        point_list = np.arange(xlist.size)
    else:
        coords = np.loadtxt(pointfile) * r_e
        point_list = np.arange(len(coords))
    if not intpol:
        fobj = pt.vlsvfile.VlsvReader(
            bulkpath + bulkprefix + ".{}.vlsv".format(str(fnr0).zfill(7))
        )
        cellids = [int(fobj.get_cellid(coord)) for coord in coords]

    data_arr = np.zeros((fnr_arr.size, point_list.size), dtype=float)

    for idx in range(fnr_arr.size):
        fnr = fnr_arr[idx]
        vlsvobj = pt.vlsvfile.VlsvReader(
            bulkpath + bulkprefix + ".{}.vlsv".format(str(fnr).zfill(7))
        )
        t_arr[idx] = vlsvobj.read_parameter("time")
        if intpol:
            data_arr[idx, :] = [
                vlsvobj.read_interpolated_variable(var, coord, operator=op)
                for coord in coords
            ]
        else:
            data_arr[idx, :] = vlsvobj.read_variable(var, operator=op, cellids=cellids)

    if filt > 0:
        data_arr = data_arr - uniform_filter1d(data_arr, size=filt, axis=0)

    XmeshXY, YmeshXY = np.meshgrid(point_list, t_arr)

    fig, ax = plt.subplots(1, 1, figsize=(8, 12), constrained_layout=True)

    im = ax.pcolormesh(
        XmeshXY,
        YmeshXY,
        data_arr,
        shading="gouraud",
        cmap=cmap,
        rasterized=True,
    )

    ax.set_xlim(point_list[0], point_list[-1])
    ax.set_ylim(t_arr[0], t_arr[-1])
    ax.set_xlabel("Point along cut", labelpad=10, fontsize=16)
    ax.set_ylabel("Time [s]", labelpad=10, fontsize=16)
    ax.set_title(var, pad=10, fontsize=16)

    cb = fig.colorbar(im, ax=ax)

    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except OSError:
            pass

    if outputdir[-1] != "/":
        outputdir += "/"

    fig.savefig(outputdir + outputname, dpi=300)
    plt.close(fig)


parser = argparse.ArgumentParser()
parser.add_argument(
    "-var", help="Variable to plot (e.g. proton/vg_rho)", type=str, required=True
)
parser.add_argument("-op", help="Operator for variable", type=str, default="pass")
parser.add_argument(
    "-fnr", nargs=2, help="First and last file number to plot", type=int, required=True
)
parser.add_argument(
    "-pointfile",
    help="A file with coordinates defining the curve to plot along",
    type=str,
)
parser.add_argument(
    "-startpoint",
    nargs=3,
    help="Start coords of straight line cut-through in Re",
    type=float,
)
parser.add_argument(
    "-endpoint",
    nargs=3,
    help="End coords of straight line cut-through in Re",
    type=float,
)
parser.add_argument(
    "-dr", help="distance between points in straight line cut-through [km]", type=float
)
parser.add_argument(
    "-bulkpath", help="Path to directory with bulk files", type=str, required=True
)
parser.add_argument(
    "-bulkprefix",
    help="Starting string of bulk file name (e.g. bulk, bulk1,bulk5)",
    type=str,
    required=True,
)
parser.add_argument("-outputdir", help="Output file directory", type=str)
parser.add_argument("-outputname", help="Output file name", type=str)
parser.add_argument(
    "-intpol",
    help="Interpolate to cut-through sample points? (True/False)",
    type=bool,
    default=False,
)
parser.add_argument(
    "-filt",
    help="Filter out slowly changing signal? (<=0: no filtering, >0: filter with specified window size)",
    type=int,
    default=-1,
)
parser.add_argument(
    "-cmap", help="Colormap to use for plot", type=str, default="viridis"
)
args = parser.parse_args()

if __name__ == "__main__":

    # if len(sys.argv) != 13:
    #     print("\n")
    #     print("This script takes 12 parameters, example usage:")
    #     print(
    #         "python cutthrough_timeseries.py var fnr0 fnr1 x0 y0 z0 x1 y1 z1 dr bulkpath bulkprefix outputname outputdir intpol filt op cmap"
    #     )
    #     print("Parameter descriptions:")
    #     print("var: Variable to plot")
    #     print("fnr0: First file number to plot")
    #     print("fnr1: Last file number to plot")
    #     # print("x0: Cut-through starting point x coordinate [Re]")
    #     # print("y0: Cut-through starting point y coordinate [Re]")
    #     # print("z0: Cut-through starting point z coordinate [Re]")
    #     # print("x1: Cut-through ending point x coordinate [Re]")
    #     # print("y1: Cut-through ending point y coordinate [Re]")
    #     # print("z1: Cut-through ending point z coordinate [Re]")
    #     # print("dr: Distance between cut-through sample points (km)")
    #     print("pointfile: File containing a list of coordinates, in Re")
    #     print("bulkpath: Path to bulk files")
    #     print("bulkprefix: Starting string of bulk file name (e.g. bulk, bulk1, bulk5)")
    #     print("outputname: Name of output file")
    #     print("outputdir: Output file directory")
    #     print("intpol: Interpolate to cut-through sample points? (True/False)")
    #     print(
    #         "filt: Filter out slowly changing signal? (<=0: no filtering, >0: filter with specified window size)"
    #     )
    #     print("op: Variable operator")
    #     print("cmap: Colormap")
    #     print("\nExiting")
    #     print("\n")
    #     raise Exception

    # (
    #     arg0,
    #     var,
    #     fnr0,
    #     fnr1,
    #     # x0,
    #     # y0,
    #     # z0,
    #     # x1,
    #     # y1,
    #     # z1,
    #     # dr,
    #     pointfile,
    #     bulkpath,
    #     bulkprefix,
    #     outputname,
    #     outputdir,
    #     intpol,
    #     filt,
    #     op,
    #     cmap,
    # ) = sys.argv

    jplots(
        var=args.var,
        fnr0=args.fnr[0],
        fnr1=args.fnr[1],
        start_coords=args.startpoint,
        end_coords=args.endpoint,
        dr=args.dr,
        bulkpath=args.bulkpath,
        bulkprefix=args.bulkprefix,
        outputdir=args.outputdir,
        outputname=args.outputname,
        intpol=args.intpol,
        filt=args.filt,
        op=args.op,
        cmap=args.cmap,
        pointfile=args.pointfile,
    )
