import sys, os
import numpy as np
from scipy.ndimage import uniform_filter1d
import pytools as pt
import matplotlib.pyplot as plt
import argparse
from pyCalculations.cutthrough import cut_through

r_e = 6.371e6

"""
    Script for plotting a cut-through timeseries, aka keogram, aka time-elongation map, etc. for a given variable, a given set of coordinates, and a range of times.

    This script takes 12 parameters, example usage:
    
    python cutthrough_timeseries.py -var <var> -fnr <fnr1> <fnr2> -bulkpath <bulkpath> -bulkprefix <bulkprefix> -outputname <outputname> -outputdir <outputdir> -filt <filt> -op <op> -cmap <cmap> -point <point1> <point2>

    
    Parameter descriptions:
        var: Variable to plot (e.g. "proton/vg_rho", "vg_b_vol")
        fnr1: First file number to plot
        fnr2: Last file number to plot
        bulkpath: Path to directory containing bulk files
        bulkprefix: Starting string of bulk file name (e.g. bulk, bulk1, bulk5)
        point1: First point on the line to plot
        point2: Last point on the line to plot
        outputname: Name of output file
        outputdir: Path to output directory
        filt: Filter out temporally slowly changing signal? (<=0: no filtering, >0: filter with specified window size), default=-1
        op: Variable operator, default="pass"
        cmap: Colormap, default="viridis"
"""


def jplots(
    var,
    fnr0,
    fnr1,
    bulkpath,
    bulkprefix,
    point1,
    point2,
    outputname,
    outputdir,
    filt=-1,
    op="pass",
    cmap="viridis",
):

    fnr_arr = np.arange(fnr0, fnr1 + 0.1, 1, dtype=int)
    t_arr = np.zeros_like(fnr_arr).astype(float)

    if bulkpath[-1] != "/":
        bulkpath += "/"

    fobj = pt.vlsvfile.VlsvReader(
        bulkpath + bulkprefix + ".{}.vlsv".format(str(fnr0).zfill(7))
    )
    cut = cut_through(fobj, point1, point2)
    cellids = cut["CellID"][:-1]
    distances = cut["distances"]
    distances = np.array(distances) / r_e

    data_arr = np.zeros((fnr_arr.size, distances.size), dtype=float)

    for idx in range(fnr_arr.size):
        fnr = fnr_arr[idx]
        vlsvobj = pt.vlsvfile.VlsvReader(
            bulkpath + bulkprefix + ".{}.vlsv".format(str(fnr).zfill(7))
        )
        t_arr[idx] = vlsvobj.read_parameter("time")
        data_arr[idx, :] = vlsvobj.read_variable(var, operator=op, cellids=cellids)
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

    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except OSError:
            pass

    if outputdir[-1] != "/":
        outputdir += "/"

    fig.savefig(outputdir + outputname, dpi=300)
    plt.close(fig)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-var", help="Variable to plot (e.g. proton/vg_rho)", type=str, required=True
    )
    parser.add_argument("-op", help="Operator for variable", type=str, default="pass")
    parser.add_argument(
        "-fnr",
        nargs=2,
        help="First and last file number to plot",
        type=int,
        required=True,
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
    parser.add_argument(
        "-point",
        nargs=2,
        help="First and last point on the line to plot",
        type=list,
        required=True,
    )
    parser.add_argument("-outputdir", help="Output file directory", type=str)
    parser.add_argument("-outputname", help="Output file name", type=str)
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

    jplots(
        var=args.var,
        fnr1=args.fnr[0],
        fnr2=args.fnr[1],
        bulkpath=args.bulkpath,
        bulkprefix=args.bulkprefix,
        point1=args.point[0],
        point2=args.point[1],
        outputdir=args.outputdir,
        outputname=args.outputname,
        filt=args.filt,
        op=args.op,
        cmap=args.cmap,
    )


if __name__ == "__main__":

    main()
