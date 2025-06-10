'''
Finds the magnetopause position by tracing streamlines of the plasma flow for two-dimensional Vlasiator runs. Needs the yt package.
'''

import numpy as np
import analysator as pt
import yt


def interpolate(streamline, x_points):
    """Interpolates a single streamline for make_magnetopause().

        :param streamline: a single streamline to be interpolated
        :param x_points: points in the x-axis to use for interpolation
        :returns: the streamline as numpy array of x,z coordinate points where the x-axis coordinates are the points given to the function
    """

    arr = np.array(streamline)

    # set arrays for interpolation
    xp = arr[:,0][::-1]
    zp = arr[:,2][::-1]

    #interpolate z coordinates
    z_points = np.interp(x_points, xp, zp, left=np.NaN, right=np.NaN)

    return np.array([x_points, z_points])


def make_streamlines(vlsvfile, streamline_seeds=None, streamline_length=40*6371000):
    """Traces streamlines of velocity field using the yt package.

        :param vlsvfile: directory and file name of .vlsv data file to use for VlsvReader
        :kword streamline_seeds: optional streamline starting points in numpy array (coordinates in meters including the y-coordinate 0.0)
        :kword streamline_length: streamline length

        :returns: streamlines as numpy array
    """
     
    # bulk file
    f = pt.vlsvfile.VlsvReader(file_name=vlsvfile)

    # get box coordinates from data
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    simext =[xmin,xmax,ymin,ymax,zmin,zmax]
    sizes = np.array([xsize,ysize,zsize])
    boxcoords=list(simext)

    cellids = f.read_variable("CellID")

    #Read the data from vlsv-file
    Vx = f.read_variable("v", operator="x")
    Vz = f.read_variable("v", operator="z")


    #Re-shape variable data
    Vxs=Vx[np.argsort(cellids)].reshape(f.get_spatial_mesh_size(), order="F")
    Vys = np.zeros_like(Vxs)
    Vzs=Vz[np.argsort(cellids)].reshape(f.get_spatial_mesh_size(), order="F")

    data=dict(Vx=Vxs,Vy=Vys,Vz=Vzs)

    #Create starting points for streamlines if they are not given
    if streamline_seeds == None:
        streamline_seeds = np.array([[20*6371000, 0 ,i] for i in  np.linspace(-5*6371000, 5*6371000, 200)])

    #streamline_seeds = np.array(streamline_seeds)
    #dataset in yt-form
    yt_dataset = yt.load_uniform_grid(
        data,
        sizes,
        bbox=np.array([[boxcoords[0], boxcoords[1]],
                    [boxcoords[2],boxcoords[3]],
                    [boxcoords[4],boxcoords[5]]]))
                                    

    #data, seeds, dictionary positions, step size
    streamlines = yt.visualization.api.Streamlines(yt_dataset, streamline_seeds,
                                                    "Vx", "Vy", "Vz", length=streamline_length, direction=1)

    #trace the streamlines with yt
    streamlines.integrate_through_volume()
    # return streamline positions
    return np.array(streamlines.streamlines)


def make_magnetopause(streamlines, end_x=-15*6371000, x_point_n=50):
    """Finds the mangetopause location based on streamlines.

        :param streams: streamlines (coordinates in m)
        :kword end_x: tail end x-coordinate (how far along the negative x-axis the magnetopause is calculated)
        :kword x_point_n: integer, how many x-axis points the magnetopause will be divided in between the subsolar point and tail
        
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
    new_streampoints = np.zeros((len(x_points), len(streamlines), 1)) # new array for keeping interpolated streamlines in form streamlines_new[x_point, streamline, z-coordinate] 

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


def find_magnetopause(vlsvfile, streamline_seeds=None, streamline_length=45*6371000, end_x=-15*6371000, x_point_n=50):
    """Finds the magnetopause position by tracing streamlines of the velocity field for 2d runs.

        :param vlsvfile: directory and file name of .vlsv data file to use for VlsvReader
        :kword streamline_seeds: optional streamline starting points in numpy array (coordinates in meters including the y-coordinate 0.0)
        :kword streamline_length: streamline length for tracing
        :kword end_x: tail end x-coordinate (how far along the negative x-axis the magnetopause is calculated)
        :kword x_point_n: integer, how many x-axis points the magnetopause will be divided in between the subsolar point and tail

        :returns:   the magnetopause position as coordinate points in numpy array
    """

    streamlines = make_streamlines(vlsvfile, streamline_seeds, streamline_length)
    magnetopause = make_magnetopause(streamlines, end_x, x_point_n)

    return magnetopause
