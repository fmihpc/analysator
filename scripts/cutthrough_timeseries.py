import sys, os
import numpy as np
from scipy.ndimage import uniform_filter1d
import pytools as pt
import matplotlib.pyplot as plt

r_e = 6.371e6


def jplots(
    var,
    fnr0,
    fnr1,
    x0,
    y0,
    z0,
    x1,
    y1,
    z1,
    dr,
    bulkpath,
    bulkprefix,
    outputname,
    outputdir="./",
    intpol=False,
    filt=-1,
    op="pass",
    cmap="viridis",
):

    dr /= r_e

    fnr_arr = np.arange(fnr0, fnr1 + 1, 0.1, dtype=int)
    t_arr = np.zeros_like(fnr_arr).astype(float)

    if bulkpath[-1] != "/":
        bulkpath += "/"

    npoints = int(np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2) / dr) + 1

    xlist = np.linspace(x0, x1, npoints, dtype=float)
    ylist = np.linspace(y0, y1, npoints, dtype=float)
    zlist = np.linspace(z0, z1, npoints, dtype=float)

    if intpol:
        coords = r_e * np.array([xlist, ylist, zlist]).T
    else:
        fobj = pt.vlsvfile.VlsvReader(
            bulkpath + bulkprefix + ".{}.vlsv".format(str(fnr0).zfill(7))
        )
        cellids = [
            int(fobj.get_cellid([xlist[idx] * r_e, ylist[idx] * r_e, 0]))
            for idx in range(xlist.size)
        ]
    point_list = np.arange(xlist.size)

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
        data_arr = uniform_filter1d(data_arr, size=filt, axis=0)

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
    ax.set_xlabel("Point along cut", labelpad=10)
    ax.set_ylabel("Time [s]", labelpad=10)

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


if __name__ == "__main__":
    (
        arg0,
        var,
        fnr0,
        fnr1,
        x0,
        y0,
        z0,
        x1,
        y1,
        z1,
        dr,
        bulkpath,
        bulkprefix,
        outputname,
        outputdir,
        intpol,
        filt,
        op,
        cmap,
    ) = sys.argv

    jplots(
        var,
        int(fnr0),
        int(fnr1),
        float(x0),
        float(y0),
        float(z0),
        float(x1),
        float(y1),
        float(z1),
        float(dr),
        bulkpath,
        bulkprefix,
        outputname,
        outputdir,
        bool(intpol),
        int(filt),
        op,
        cmap,
    )
