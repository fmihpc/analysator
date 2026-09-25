'''
Finds the magnetopause position by tracing streamlines of the plasma flow for two-dimensional Vlasiator runs.

Streamlines are integrated with scipy using RegularGridInterpolator for the velocity field and solve_ivp for the ODE integration
'''

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.integrate import solve_ivp
import analysator as pt


def interpolate(streamline, x_points):
    """Interpolates a single streamline for make_magnetopause().

        :param streamline: a single streamline to be interpolated
        :param x_points: points in the x-axis to use for interpolation
        :returns: the streamline as numpy array of x,z coordinate points where the x-axis coordinates are the points given to the function
    """

    arr = np.array(streamline)

    # Drop points beyond which the streamline has no data
    valid = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 2])
    arr = arr[valid]

    x_points = np.asarray(x_points, dtype=float)
    if arr.shape[0] == 0:
        return np.array([x_points, np.full_like(x_points, np.nan)])

    # set arrays for interpolation
    xp = arr[:,0][::-1]
    zp = arr[:,2][::-1]

    # interpolate z coordinates
    z_points = np.interp(x_points, xp, zp, left=np.nan, right=np.nan)

    return np.array([x_points, z_points])


def trace_streamline(x0, z0, vx_interp, vz_interp, bounds, length, max_step, n_points, rtol, atol, min_speed):
    """Integrates one single streamline forward starting at (x0, z0).

        :param x0, z0: starting point, m
        :param vx_interp, vz_interp: RegularGridInterpolator for the two velocity components (must return nan outside their grid)
        :param bounds: (xmin, xmax, zmin, zmax) -- valid range of the interpolators, m
        :param length: maximum arc length to integrate, m
        :param max_step: maximum solve_ivp step size, m
        :param n_points: number of points to sample along the streamline
        :param rtol, atol: solve_ivp error tolerances
        :param min_speed: velocity magnitude (m/s) below which the local flow direction is treated as undefined -- guards against 0/0 in regions of exactly zero velocity (e.g. an inner boundary)

        :returns: (n_points, 3) array of x, y(=0), z points. Points beyond where the streamline left the domain are nan.
    """
    xmin, xmax, zmin, zmax = bounds

    def rhs(s, state):
        x, z = state
        vx = float(vx_interp((x, z)))
        vz = float(vz_interp((x, z)))
        speed = np.hypot(vx, vz)
        if not np.isfinite(speed) or speed < min_speed:
            return (0.0, 0.0)
        return (vx / speed, vz / speed)

    def exit_domain(s, state):
        x, z = state
        return min(x-xmin, xmax-x, z-zmin, zmax-z)
    exit_domain.terminal = True
    exit_domain.direction = -1

    sol = solve_ivp(rhs, [0.0, length], [x0, z0], events=exit_domain, dense_output=True, max_step=max_step, rtol=rtol, atol=atol)

    s_eval = np.linspace(0.0, length, n_points)
    s_max = sol.t[-1]

    out = np.full((n_points, 3), np.nan)
    valid = s_eval <= s_max
    pts = sol.sol(s_eval[valid])
    out[valid, 0] = pts[0]
    out[valid, 1] = 0.0
    out[valid, 2] = pts[1]
    return out


def make_streamlines(vlsvfile, streamline_seeds=None, seeds_n=200, seeds_x0=20*6371000, seeds_range=(-5*6371000, 5*6371000), streamline_length=40*6371000, *, max_step=None, n_points=400, rtol=1e-6, atol=1.0, min_speed=1e-3):
    """Traces streamlines of the velocity field.

        :param vlsvfile: directory and file name of .vlsv data file to use for VlsvReader
        :kwarg streamline_seeds: optional streamline starting points in numpy array (coordinates in meters including the y-coordinate 0.0)
        :kwarg seeds_n: instead of streamline_seeds provide a number of streamlines to be traced
        :kwarg seeds_x0: instead of streamline_seeds provide an x-coordinate for streamline starting points
        :kwarg seeds_range: instead of streamline_seeds provide [min, max] range to use for streamline starting point z-coordinates
        :kwarg streamline_length: streamline length
        :kwarg max_step: maximum integration step, m (default: one grid cell of the run)
        :kwarg n_points: number of points sampled along each returned streamline
        :kwarg rtol, atol: scipy.integrate.solve_ivp error tolerances
        :kwarg min_speed: velocity magnitude (m/s) treated as numerically zero

        :returns: streamlines as numpy array, shape (n_seeds, n_points, 3). Points beyond which a streamline left the simulation domain are nan.
    """

    # bulk file
    f = pt.vlsvfile.VlsvReader(file_name=vlsvfile)

    # get box coordinates from data
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    mesh_size = f.get_spatial_mesh_size()
    [xsize, ysize, zsize] = mesh_size

    cellids = f.read_variable("CellID")

    # Read the data from vlsv-file
    Vx = f.read_variable("v", operator="x")
    Vz = f.read_variable("v", operator="z")

    # Re-shape variable data
    order = np.argsort(cellids)
    Vxs = Vx[order].reshape(mesh_size, order="F")
    Vzs = Vz[order].reshape(mesh_size, order="F")

    # this routine is for 2D (x-z plane) runs: a single cell in y
    Vxs2d = Vxs[:, 0, :]
    Vzs2d = Vzs[:, 0, :]

    # cell-centre coordinates of the (uniform) spatial mesh
    dx = (xmax - xmin) / xsize
    dz = (zmax - zmin) / zsize
    x_coords = xmin + dx * (np.arange(xsize) + 0.5)
    z_coords = zmin + dz * (np.arange(zsize) + 0.5)

    vx_interp = RegularGridInterpolator((x_coords, z_coords), Vxs2d, method="linear", bounds_error=False, fill_value=np.nan)
    vz_interp = RegularGridInterpolator((x_coords, z_coords), Vzs2d, method="linear", bounds_error=False, fill_value=np.nan)
    bounds = (x_coords[0], x_coords[-1], z_coords[0], z_coords[-1])

    if max_step is None:
        max_step = min(dx, dz)

    # Create starting points for streamlines if they are not given
    if streamline_seeds is None:
        streamline_seeds = np.array([[seeds_x0, 0.0, i] for i in np.linspace(seeds_range[0], seeds_range[1], seeds_n)])
    else:
        streamline_seeds = np.asarray(streamline_seeds)

    streamlines = np.array([
        trace_streamline(x0, z0, vx_interp, vz_interp, bounds, streamline_length, max_step, n_points, rtol, atol, min_speed)
        for x0, _, z0 in streamline_seeds
    ])

    return streamlines


def make_magnetopause(streamlines, end_x=-15*6371000, x_point_n=50):
    """Finds the mangetopause location based on streamlines.

        :param streams: streamlines (coordinates in m)
        :kwarg end_x: tail end x-coordinate (how far along the negative x-axis the magnetopause is calculated)
        :kwarg x_point_n: integer, how many x-axis points the magnetopause will be divided in between the subsolar point and tail

        :returns:   the magnetopause position as coordinate points in numpy array
    """

    RE = 6371000

    streampoints = np.reshape(streamlines, (streamlines.shape[0]*streamlines.shape[1], 3)) #all the points in one array

    ## find the subsolar dayside point in the positive x-axis
    ## do this by finding a stremline point on positive x axis closest to the Earth
    x_axis_points = streampoints[(streampoints[:,2]< RE) & (streampoints[:,2]> -RE) & (streampoints[:,0]> 0)]
    subsolar_x =np.min(x_axis_points[:,0])

    ## define points in the x axis where to find magnetopause points on the yz-plane
    x_points = np.linspace(subsolar_x, end_x, x_point_n)

    ## interpolate more exact points for streamlines at exery x_point
    # new array for keeping interpolated streamlines in form streamlines_new[x_point, streamline, z-coordinate]
    new_streampoints = np.zeros((len(x_points), len(streamlines), 1))

    for i,stream in enumerate(streamlines):
        interpolated_streamline = interpolate(stream, x_points)
        for j in range(0, len(x_points)):
            new_streampoints[j, i,:] = interpolated_streamline[1,j]


    ## start making the magnetopause
    ## in each x_point, find the closest streamline to x-axis in the positive and negative z-axis

    pos_z_mpause = np.zeros((len(x_points), 2))
    neg_z_mpause =  np.zeros((len(x_points), 2))

    for i, x_point in enumerate(x_points):
        pos = new_streampoints[i, new_streampoints[i,:] > 0]
        neg = new_streampoints[i, new_streampoints[i,:] < 0]

        if (pos.size == 0) or (neg.size == 0):
            raise ValueError('No streamlines found for x axis point, try adding streamlines or checking the x_points')

        # find points closest to x-axis and save found points
        pos_z_mpause[i] = [x_point, pos[pos.argmin()]]
        neg_z_mpause[i] = [x_point, neg[neg.argmax()]]

    magnetopause = np.concatenate((pos_z_mpause[::-1], np.array([[subsolar_x, 0]]),  neg_z_mpause))

    return magnetopause


def find_magnetopause_sw_streamline_2d(vlsvfile, streamline_seeds=None, seeds_n=200, seeds_x0=20*6371000, seeds_range=(-5*6371000, 5*6371000), streamline_length=45*6371000, end_x=-15*6371000, x_point_n=50, *, max_step=None, n_points=400, rtol=1e-6, atol=1.0, min_speed=1e-3):
    """Finds the magnetopause position by tracing streamlines of the velocity field for 2d runs.

        :param vlsvfile: directory and file name of .vlsv data file to use for VlsvReader
        :kwarg streamline_seeds: optional streamline starting points in numpy array (coordinates in meters including the y-coordinate 0.0)
        :kwarg seeds_n: instead of streamline_seeds provide a number of streamlines to be traced
        :kwarg seeds_x0: instead of streamline_seeds provide an x-coordinate for streamline starting points
        :kwarg seeds_range: instead of streamline_seeds provide [min, max] range to use for streamline starting point z-coordinates
        :kwarg streamline_length: streamline length for tracing
        :kwarg end_x: tail end x-coordinate (how far along the negative x-axis the magnetopause is calculated)
        :kwarg x_point_n: integer, how many x-axis points the magnetopause will be divided in between the subsolar point and tail
        :kwarg max_step, n_points, rtol, atol, min_speed: passed through to make_streamlines(); see there

        :returns:   the magnetopause position as coordinate points in numpy array
    """

    streamlines = make_streamlines(vlsvfile, streamline_seeds, seeds_n, seeds_x0, seeds_range, streamline_length, max_step=max_step, n_points=n_points, rtol=rtol, atol=atol, min_speed=min_speed)
    magnetopause = make_magnetopause(streamlines, end_x, x_point_n)

    return magnetopause
