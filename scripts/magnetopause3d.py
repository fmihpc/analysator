'''
Finds the magnetopause position by tracing streamlines of the plasma flow for three-dimensional Vlasiator runs.
'''

import matplotlib.pyplot as plt
import numpy as np
import analysator as pt

def cartesian_to_polar(cartesian_coords): # for segments of plane
    """Converts cartesian coordinates to polar (for the segments of the yz-planes).

        :param cartesian_coords: (y,z) coordinates as list or array
        :returns: the polar coordinates r, phi (angle in degrees)
    """
    y,z = cartesian_coords[0], cartesian_coords[1]
    r = np.sqrt(z**2 + y**2)
    phi = np.arctan2(z, y)
    phi = np.rad2deg(phi) #angle in degrees
    if phi < 0: phi = phi+360
    return(r, phi)

def polar_to_cartesian(r, phi):
    """Converts polar coordinates of the yz-plane to cartesian coordinates.

        :param r: radius of the segment (distance from the x-axis)
        :param phi: the angle coordinate in degrees
        :returns:  y, z -coordinates in cartesian system
    """
    phi = np.deg2rad(phi)
    y = r * np.cos(phi)
    z = r * np.sin(phi)
    return(y, z)


def interpolate(streamline, x_points):
    """Interpolates a single streamline for make_magnetopause(). 

        :param streamline: a single streamline to be interpolated
        :param x_points: points in the x-axis to use for interpolation
        :returns: the streamline as numpy array of coordinate points where the x-axis coordinates are the points given to the function
    """
    arr = np.array(streamline)

    # set arrays for interpolation
    xp = arr[:,0][::-1]
    yp = arr[:,1][::-1]
    zp = arr[:,2][::-1]

    #interpolate z from xz-plane
    z_points = np.interp(x_points, xp, zp, left=np.NaN, right=np.NaN)
    #interpolate y from xy-plane
    y_points = np.interp(x_points, xp, yp, left=np.NaN, right=np.NaN)


    return np.array([x_points, y_points, z_points])




def make_surface(coords):
    '''Defines a surface constructed of input coordinates as triangles.

        :param coords: points that make the surface
        :returns: list of verts and vert indices of surface triangles as numpy arrays

        input coordinates must be in form [...[c11, c21, c31, ... cn1],[c12, c22, c32, ... cn2],...
        where cij = [xij, yij, zij], i marks sector, j marks yz-plane (x_coord) index 

        Three points make a triangle, triangles make the surface.
        For every two planes next to each other:
        - take every other point from plane1, every other from plane2 (in order!)
        - from list of points: every three points closest to each other make a surface

        Example:
        plane 1: [v1, v2, v3, v4]
        plane 2: [v5, v6, v7, v8]

        -> list: [v1, v5, v2, v6, v3,...]
        -> triangles:
        v1 v5 v2
        v5 v2 v6
        v2 v6 v3
        .
        .
        .

    '''
    verts = [] #points
    faces = [] #surface triangles

    slices_in_plane = len(coords[0])
    planes = len(coords)

    #get points
    for plane in coords:
        for vert in plane:
            verts.append(vert)

    #get triangle (face) indices
    #Let's consider the area between the first two planes
    #ind1, ind2, ind3 for triangle indices
    ind1 = 0 #starts from plane 1
    ind2 = slices_in_plane #starts from plane 2
    ind3 = 1 #starts from plane 1
    first_triangles = []

    while len(first_triangles) < (2*slices_in_plane):
        first_triangles.append([ind1, ind2, ind3])
        ind1 = ind2
        ind2 = ind3
        if (ind3 == (slices_in_plane*2)-1): #max index, go over to plane 1 first index
            ind3 = 0
        elif (ind3 == 0): #last round, go to plane 2 first index
            ind3 = slices_in_plane
        else:
            ind3 = ind1 + 1

    first_triangles = np.array(first_triangles)
    
    #Now the rest of the triangles are just the first triangles + (index of area * slices_in_plane)
    for area_index in range(planes-1):
        next_triangles = [x + slices_in_plane*area_index for x in first_triangles]
        faces.extend(next_triangles)

    return np.array(verts), np.array(faces)


def make_streamlines(vlsvfile, streamline_seeds=None, dl=2e6, iterations=200):
    """Traces streamlines of velocity field from outside the magnetosphere to magnetotail.

        :param vlsvfile: directory and file name of .vlsv data file to use for VlsvReader

        :kword streamline_seeds: optional streamline starting points in numpy array
        :kword dl: streamline iteration step length in m
        :kword iterations: int, number of iteration steps

        :returns: streamlines as numpy array
    """

    f = pt.vlsvfile.VlsvReader(file_name=vlsvfile)

    # Create streamline starting points if they have not been given
    if streamline_seeds == None:
        streamline_seeds = np.zeros([25**2, 3])

        t = np.linspace(-5*6371000, 5*6371000, 25)
        k = 0
        for i in t:
            for j in t:
                streamline_seeds[k] = [20*6371000, i, j]
                k = k+1

    # Trace the streamlines
    streams = pt.calculations.static_field_tracer_3d(
        vlsvReader=f,
        seed_coords=streamline_seeds,
        max_iterations=iterations,
        dx=dl,
        direction='+',
        grid_var='vg_v')

    return streams


def make_magnetopause(streams, end_x=-15*6371000, x_point_n=50, sector_n=36):
    """Finds the mangetopause location based on streamlines.

        :param streams: streamlines (coordinates in m)

        :kword end_x: tail end x-coordinate (how far along the negative x-axis the magnetopause is calculated)
        :kword x_point_n: integer, how many x-axis points the magnetopause will be divided in between the subsolar point and tail
        :kword sector_n: integer, how many sectors the magnetopause will be divided in on each yz-plane

        :returns:   the magnetopause position as coordinate points in numpy array, form [...[c11, c21, c31, ... cn1],[c12, c22, c32, ... cn2],...
                    where cij = [xij, yij, zij], i marks sector, j marks yz-plane (x-coordinate) index 
    """

    RE = 6371000

    #streams = streams*(1/RE) # streamlines in rE
    streampoints = np.reshape(streams, (streams.shape[0]*streams.shape[1], 3)) #all the points in one array)
    
    ## find the subsolar dayside point in the x-axis
    ## do this by finding a streamline point on positive x axis closest to the Earth
    # streampoints closer than ~1 rE to positive x-axis:
    x_axis_points = streampoints[(streampoints[:,1]<RE) & (streampoints[:,2]<RE) & (streampoints[:,1]>-RE) & (streampoints[:,2]>-RE) & (streampoints[:,0]>0) & (streampoints[:,0]>0)] 
    subsolar_x =np.min(x_axis_points[:,0])


    ## define points in the x axis where to find magnetopause points on the yz-plane
    x_points = np.linspace(subsolar_x, end_x, x_point_n)
    
    ## interpolate more exact points for streamlines at exery x_point
    new_streampoints = np.zeros((len(x_points), len(streams), 2)) # new array for keeping interpolated streamlines in form new_streampoints[x_point, streamline, y and z -coordinates] 
    i=0 # streamline
    for stream in streams:
        interpolated_streamline = interpolate(stream, x_points)

        if type(interpolated_streamline) is np.ndarray: # don't use 'discarded' streamlines, see function interpolate()
            for j in range(0, len(x_points)):
                y,z = interpolated_streamline[1,j], interpolated_streamline[2,j]
                new_streampoints[j, i,:] = np.array([y,z])
        i += 1



    ## create a list of streamline points in polar coordinates (for every x_point)
    polar_coords = np.zeros_like(new_streampoints)
    for i in range(0,new_streampoints.shape[0]):
        for j in range(0,new_streampoints.shape[1]):
            polar_coords[i,j,:] = cartesian_to_polar(new_streampoints[i,j])


    ## now start making the magnetopause
    ## in each x_point, divide the plane into sectors and look for the closest streamline to x-axis in the sector

    ## if given sector number isn't divisible by 4, make it so because we want to have magnetopause points at exactly y=0 and z=0 for 2d slices of the whole thing
    while sector_n%4 != 0:
        sector_n +=1

    sector_width = 360/sector_n
    magnetopause = np.zeros((len(x_points), sector_n, 3))

    for i,x_point in enumerate(x_points): #loop over every chosen x-axis point
        # divide the yz-plane into sectors
        for j, mean_sector_angle in enumerate(np.arange(0, 360, sector_width)):
            min_angle = mean_sector_angle-sector_width/2
            max_angle = mean_sector_angle+sector_width/2

            # find points that are in the sector
            if mean_sector_angle == 0: # special case as the first sector needs streamlines around phi=0
                min_angle = min_angle+360
                # divide into phi<360 and phi>0
                sector1 = polar_coords[i, (polar_coords[i,:,1] <= 360)*(polar_coords[i,:,1] > min_angle)]
                sector2 = polar_coords[i, (polar_coords[i,:,1] <= max_angle)*(polar_coords[i,:,1] >= 0)]
                sector_points = np.concatenate((sector1, sector2))

            else:
                sector_points = polar_coords[i, (polar_coords[i,:,1] <= max_angle)*(polar_coords[i,:,1] > min_angle)]

            # discard 'points' with r=0 and check that there's at least one streamline point in the sector
            sector_points = sector_points[sector_points[:,0] != 0.0]
            if sector_points.size == 0:
                raise ValueError('No streamlines found in the sector')

            # find the points closest to the x-axis
            closest_point_radius = sector_points[sector_points[:,0].argmin(), 0] # smallest radius
            
            # return to cartesian coordinates and save as a magnetopause point at the middle of the sector
            y,z = polar_to_cartesian(closest_point_radius, mean_sector_angle)
            magnetopause[i,j,:] = [x_point, y, z]


    # make a tip point for the magnetopause for prettier 3d plots 
    tip = np.array([subsolar_x, 0, 0])
    tips = np.tile(tip, (magnetopause.shape[1],1))
    magnetopause = np.vstack(([tips], magnetopause))
    
    return magnetopause


def find_magnetopause(vlsvfile, streamline_seeds=None, dl=2e6, iterations=200, end_x=-15*6371000, x_point_n=50, sector_n=36):
    """Finds the magnetopause position by tracing streamlines of the velocity field.

        :param vlsvfile: path to .vlsv bulk file to use for VlsvReader
        :kword streamline_seeds: optional streamline starting points in numpy array
        :kword dl: streamline iteration step length in m
        :kword iterations: int, number of iteration steps
        :kword end_x: tail end x-coordinate (how far along the x-axis the magnetopause is calculated)
        :kword x_point_n: integer, how many x-axis points the magnetopause will be divided in between the subsolar point and tail
        :kword sector_n: integer, how many sectors the magnetopause will be divided in on each yz-plane

        :returns:   vertices, faces of the magnetopause triangle mesh as numpy arrays
    """

    streams = make_streamlines(vlsvfile, streamline_seeds, dl, iterations)
    magnetopause = make_magnetopause(streams, end_x, x_point_n, sector_n)
    vertices, faces = make_surface(magnetopause)

    return vertices, faces
