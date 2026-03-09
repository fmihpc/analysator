import sys, os
import numpy as np
from scipy.ndimage import uniform_filter1d
import analysator as pt
import matplotlib.pyplot as plt
import argparse
from analysator.calculations.cutthrough import cut_through
from analysator.calculations.lineout import lineout

r_e = 6.371e6

"""
    Script for plotting a cut-through timeseries, aka keogram, aka time-elongation map, etc. for a given variable, a given set of coordinates, and a range of times.

    This script takes 12 parameters, example usage:
    
    python cutthrough_timeseries.py -var <var> -fnr <fnr1> <fnr2> -bulkpath <bulkpath> -bulkprefix <bulkprefix> -point <point1_x> <point1_y> <point1_z> <point2_x> <point2_y> <point2_z> -outputname <outputname> -outputdir <outputdir> -filt <filt> -op <op> -cmap <cmap> -re <re> -npoints <npoints> -interpolation_order <interpolation_order

    
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
        re: Input points given in Earth radii? default=False
        npoints: Number of points in line, default=100
        interpolation_order: Order of interpolation (0 or 1), default=1
"""


def jplots(
    var,
    fnr1,
    fnr2,
    bulkpath,
    bulkprefix,
    point1,
    point2,
    outputname,
    outputdir,
    filt=-1,
    op="pass",
    cmap="viridis",
    re=False,
    npoints=100,
    interpolation_order=1,
    draw=True,
):

    fnr_arr = np.arange(fnr1, fnr2 + 0.1, 1, dtype=int)
    t_arr = np.zeros_like(fnr_arr).astype(float)

    fobj = pt.vlsvfile.VlsvReader(
        os.path.join(bulkpath , bulkprefix , ".{}.vlsv".format(str(fnr1).zfill(7)))
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
            os.path.join(bulkpath , bulkprefix , ".{}.vlsv".format(str(fnr).zfill(7)))
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
        if not os.path.exists(outputdir):
            try:
                os.makedirs(outputdir)
            except OSError:
                pass

        fig.savefig(os.path.join(outputdir,outputname), dpi=300)
        plt.close(fig)
    else:
        return (fig,ax,XmeshXY,YmeshXY,data_arr)


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
        nargs=6,
        help="x y z of first and x y z of last point on the line to plot",
        type=float,
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
    parser.add_argument(
        "-re", help="Are input points in Earth radii?", type=bool, default=False
    )
    parser.add_argument(
        "-npoints", help="Number of points in line", type=int, default=100
    )
    parser.add_argument(
        "-interpolation_order",
        help="Order of interpolation (0 or 1)",
        type=int,
        default=1,
    )
    args = parser.parse_args()

    jplots(
        var=args.var,
        fnr1=args.fnr[0],
        fnr2=args.fnr[1],
        bulkpath=args.bulkpath,
        bulkprefix=args.bulkprefix,
        point1=[args.point[0], args.point[1], args.point[2]],
        point2=[args.point[3], args.point[4], args.point[5]],
        outputdir=args.outputdir,
        outputname=args.outputname,
        filt=args.filt,
        op=args.op,
        cmap=args.cmap,
        re=args.re,
        npoints=args.npoints,
        interpolation_order=args.interpolation_order,
        draw=False
    )


if __name__ == "__main__":

    main()
