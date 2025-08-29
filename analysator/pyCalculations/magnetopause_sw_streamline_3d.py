'''
Finds the magnetopause position by tracing streamlines of the plasma flow for three-dimensional Vlasiator runs.
'''

import matplotlib.pyplot as plt
import numpy as np
import analysator as pt
import vtk

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


def cartesian_to_spherical(cartesian_coords):
    """ Cartesian to spherical coordinates (Note: axes are swapped!)"""
    x, y, z = cartesian_coords[0], cartesian_coords[1], cartesian_coords[2]
    r = np.sqrt(x**2+y**2+z**2)
    theta = np.arctan2(np.sqrt(y**2+z**2), x) #inclination
    phi = np.arctan2(z, y)+np.pi # azimuth

    return [r, theta, phi]



def spherical_to_cartesian(sphe_coords):
    r, theta, phi = sphe_coords[0], sphe_coords[1], sphe_coords[2]
    y = r*np.sin(theta)*np.cos(phi)
    z = r*np.sin(theta)*np.sin(phi)
    x = r*np.cos(theta)
    return [x, y, z]


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
    z_points = np.interp(x_points, xp, zp, left=np.nan, right=np.nan)
    #interpolate y from xy-plane
    y_points = np.interp(x_points, xp, yp, left=np.nan, right=np.nan)


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

    #### test area ####
    # From last triangles remove every other triangle
    # (a single subsolar point -> last triangles are actual triangles instead of rectangles sliced in two)
    #removed = 0
    #for i in range(len(faces)-slices_in_plane*2, len(faces)):
    #    if i%2!=0:
    #        faces.pop(i-removed)
    #        removed += 1

    # From last triangles remove every other triangle
    # (a single subsolar point -> last triangles are actual triangles instead of rectangles sliced in two)
    # Also fix the last triangles so that they only point to one subsolar point and have normals towards outside
    #subsolar_index = int(len(verts)-slices_in_plane)

    #for i,triangle in enumerate(reversed(faces)):
    #    if i > (slices_in_plane): # faces not in last plane (we're going backwards) 
    #        break

    #    faces[len(faces)-i-1] = np.clip(triangle, a_min=0, a_max=subsolar_index)

    # this would remove duplicate subsolar points from vertices but makes 2d slicing harder
    #verts = verts[:int(len(verts)-slices_in_plane+1)]

    # Change every other face triangle (except for last slice triangles) normal direction so that all face the same way (hopefully)
    #faces = np.array(faces)
    #for i in range(len(faces)-int(slices_in_plane)):
    #    if i%2!=0:
    #        faces[i,1], faces[i,2] =  faces[i,2], faces[i,1]

    ###################

    # Change every other face triangle normal direction so that all face the same way
    faces = np.array(faces)
    for i in range(len(faces)):
        if i%2==0:
            faces[i,1], faces[i,2] =  faces[i,2], faces[i,1]

    return np.array(verts), faces


def make_vtk_surface(vertices, faces):
    """Makes a vtk DataSetSurfaceFilter from vertices and faces of a triangulated surface.

        :param vertices: vertex points of triangles in shape [[x0, y0, z0], [x1, y1, z1],...]
        :param faces: face connectivities as indices of vertices so that triangle normals point outwards
        :returns: a vtkDataSetSurfaceFilter object
    """

    points = vtk.vtkPoints()

    for vert in vertices:
        points.InsertNextPoint(vert)

    # make vtk PolyData object

    triangles = vtk.vtkCellArray()
    for face in faces:
        triangle = vtk.vtkTriangle()
        triangle.GetPointIds().SetId(0, face[0])
        triangle.GetPointIds().SetId(1, face[1])
        triangle.GetPointIds().SetId(2, face[2])
        triangles.InsertNextCell(triangle)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.Modified()
    polydata.SetPolys(triangles)
    polydata.Modified()

    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(polydata)
    surface.Update()

    return surface

def streamline_stopping_condition(vlsvReader, points, value):
   [xmin, ymin, zmin, xmax, ymax, zmax] = vlsvReader.get_spatial_mesh_extent()
   x = points[:, 0]
   y = points[:, 1]
   z = points[:, 2]
   beta_star = vlsvReader.read_interpolated_variable("vg_beta_star", points)
   return (x < xmin)|(x > xmax) | (y < ymin)|(y > ymax) | (z < zmin)|(z > zmax)|(value[:,0] > 0)| (beta_star < 0.4)


def make_streamlines(vlsvfile, streamline_seeds=None, seeds_n=25, seeds_x0=20*6371000, seeds_range=[-5*6371000, 5*6371000],  dl=2e6, iterations=200):
    """Traces streamlines of velocity field from outside the magnetosphere to magnetotail.
        Stopping condition for when streamlines turn sunwards, go out of box, or hit beta* below 0.4 region

        :param vlsvfile: directory and file name of .vlsv data file to use for VlsvReader

        :kword streamline_seeds: optional streamline starting points in numpy array
        :kword seeds_n: instead of streamline_seeds provide a number of streamlines to be traced 
        :kword seeds_x0: instead of streamline_seeds provide an x-coordinate for streamline starting points 
        :kword seeds_range: instead of streamline_seeds provide [min, max] range to use for streamline starting point coordinates (both y- and z-directions use the same range)
        :kword dl: streamline iteration step length in m
        :kword iterations: int, number of iteration steps

        :returns: streamlines as numpy array
    """

    f = pt.vlsvfile.VlsvReader(file_name=vlsvfile)

    # Create streamline starting points if needed
    if streamline_seeds == None:
        streamline_seeds = np.zeros([seeds_n**2, 3])

        t = np.linspace(seeds_range[0], seeds_range[1], seeds_n)
        k = 0
        for i in t:
            for j in t:
                streamline_seeds[k] = [seeds_x0, i, j]
                k = k+1

    # Trace the streamlines
    streams = pt.calculations.static_field_tracer_3d(
        vlsvReader=f,
        seed_coords=streamline_seeds,
        max_iterations=iterations,
        dx=dl,
        direction='+',
        grid_var='vg_v',
        stop_condition=streamline_stopping_condition
        )

    return streams


def make_magnetopause(streams, end_x=-15*6371000, x_point_n=50, sector_n=36, ignore=0):
    """Finds the magnetopause location based on streamlines.

        :param streams: streamlines (coordinates in m)

        :kword end_x: tail end x-coordinate (how far along the negative x-axis the magnetopause is calculated)
        :kword x_point_n: integer, how many x-axis points the magnetopause will be divided in between the subsolar point and tail
        :kword sector_n: integer, how many sectors the magnetopause will be divided in on each yz-plane

        :returns:   the magnetopause position as coordinate points in numpy array, form [...[c11, c21, c31, ... cn1],[c12, c22, c32, ... cn2],...
                    where cij = [xij, yij, zij], i marks sector, j marks yz-plane (x-coordinate) index 
    """

    RE = 6371000
    
    ## if given sector number isn't divisible by 4, make it so because we want to have magnetopause points at exactly y=0 and z=0 for 2d slices of the whole thing
    while sector_n%4 != 0:
        sector_n +=1

    #streams = streams*(1/RE) # streamlines in rE
    streampoints = np.reshape(streams, (streams.shape[0]*streams.shape[1], 3)) #all the points in one array)
    
    ## find the subsolar dayside point in the x-axis
    ## do this by finding a streamline point on positive x axis closest to the Earth
    # streampoints closer than ~1 rE to positive x-axis:
    x_axis_points = streampoints[(streampoints[:,1]<RE) & (streampoints[:,2]<RE) & (streampoints[:,1]>-RE) & (streampoints[:,2]>-RE) & (streampoints[:,0]>0) & (streampoints[:,0]>0)] 
    if ignore == 0:
        subsolar_x =np.min(x_axis_points[:,0])
    else:
        subsolar_x = np.partition(x_axis_points[:,0], ignore)[ignore] # take the nth point as subsolar point

    # divide the x point numbers between x > 0 (radial) an x < 0 (yz-planes) by ratio
    dayside_x_point_n =  int((subsolar_x/np.abs(end_x))*x_point_n)

    ### dayside magnetopause ###
    # for x > 0, look for magnetopause radially
    dayside_points = streampoints[streampoints[:,0] > 0]

    phi_slices = sector_n
    theta_slices = dayside_x_point_n
    phi_step = 2*np.pi/phi_slices
    theta_step = np.pi/(2*theta_slices) # positive x-axis only

    def cartesian_to_spherical_grid(cartesian_coords):
        # Only for x > 0
        r, theta, phi = cartesian_to_spherical(cartesian_coords)
        theta_idx = int(theta/theta_step)
        #phi_idx = int(phi/phi_step) # off by half grid cell in relation to x<0 grid
        if (phi <phi_step*0.5) and (phi > (0-phi_step*0.5)):
            phi_idx = 0
        else:
            phi_idx = round(phi/phi_step)

        return [theta_idx, phi_idx, r]

    def grid_mid_point(theta_idx, phi_idx):
        return (theta_idx+0.5)*theta_step, (phi_idx)*phi_step
    
    # make a dictionary based on spherical areas by theta index and phi index
    sph_points = {}
    for point in dayside_points:
        sph_gridpoint = cartesian_to_spherical_grid(point)
        idxs = (sph_gridpoint[0],sph_gridpoint[1]) # key is (theta_i, phi_i)
        if idxs not in sph_points:
            sph_points[idxs] = [sph_gridpoint[2]]
        else:
            sph_points[idxs].append(sph_gridpoint[2])
    
    # dayside magetopause from subsolar point towards origo
    dayside_magnetopause = np.zeros((theta_slices, phi_slices, 3))
    for ring_idx in range(theta_slices):
        ring_points = np.zeros((phi_slices, 3))
        for phi_idx in range(phi_slices):
            if (ring_idx, phi_idx) in sph_points:
                if ignore == 0:
                    nth_min_r = np.min(np.array(sph_points[(ring_idx, phi_idx)])) # point in area with smallest radius
                else:
                    nth_min_r = np.partition(np.array(sph_points[(ring_idx, phi_idx)]), ignore)[ignore]
            else:
                print("no ", (ring_idx, phi_idx))
                exit() # something wrong
            midpoint_theta, midpoint_phi = grid_mid_point(ring_idx, phi_idx) # the smallest radius will be assigned to a point in the middle-ish of the area
            ring_points[phi_idx] = np.array(spherical_to_cartesian([nth_min_r, midpoint_theta, midpoint_phi]))
            
        dayside_magnetopause[ring_idx] = ring_points


    ### x < 0 magnetopause ###
    # rest: look for magnetopause in yz-planes
    ## define points in the x axis where to find magnetopause points on the yz-plane
    x_points = np.linspace(0.0, end_x, x_point_n-dayside_x_point_n)
    
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
            if ignore == 0:
                nth_closest_point_radius = sector_points[sector_points[:,0].argmin(), 0] # smallest radius
            else:
                nth_closest_point_radius = np.partition(sector_points[:,0], ignore)[ignore]

            #closest_point_radius = np.median(sector_points[:,0]) # median radius
            #if x_point < subsolar_x-2e6:
            #        closest_point_radius = np.median(sector_points[:,0]) # median radius
            #else:
            #    closest_point_radius = sector_points[sector_points[:,0].argmin(), 0] # smallest radius
        
            # return to cartesian coordinates and save as a magnetopause point at the middle of the sector
            y,z = polar_to_cartesian(nth_closest_point_radius, mean_sector_angle)
            magnetopause[i,j,:] = [x_point, y, z]


    # make a tip point for the magnetopause for prettier 3d plots 
    tip = np.array([subsolar_x, 0, 0])
    tips = np.tile(tip, (magnetopause.shape[1],1))
    magnetopause = np.vstack(([tips], dayside_magnetopause, magnetopause))
    
    return magnetopause


def find_magnetopause_sw_streamline_3d(vlsvfile, streamline_seeds=None, seeds_n=25, seeds_x0=20*6371000, seeds_range=[-5*6371000, 5*6371000], dl=2e6, iterations=200, end_x=-15*6371000, x_point_n=50, sector_n=36, ignore=0):
    """Finds the magnetopause position by tracing streamlines of the velocity field.

        Note: there may be a slight jump at x=0. This may be due to difference in methods (pointcloud vs. interpolation).
            If the subsolar point area looks off then more streamlines near the x-axis are needed.


        :param vlsvfile: path to .vlsv bulk file to use for VlsvReader
        :kword streamline_seeds: optional streamline starting points in numpy array
        :kword seeds_n: instead of streamline_seeds provide a number of streamlines to be traced 
        :kword seeds_x0: instead of streamline_seeds provide an x-coordinate for streamline starting points 
        :kword seeds_range: instead of streamline_seeds provide [min, max] range to use for streamline starting point coordinates (both y- and z-directions use the same range)
        :kword dl: streamline iteration step length in m
        :kword iterations: int, number of iteration steps
        :kword end_x: tail end x-coordinate (how far along the x-axis the magnetopause is calculated)
        :kword x_point_n: integer, how many parts the magnetopause will be divided in between the subsolar point and tail end
        :kword sector_n: integer, how many sectors the magnetopause will be divided in on each yz-plane/radial sector
        :kword ignore: how many inner streamlines will be ignored when calculating the magnetopause

        :returns:   vertices, surface where vertices are numpy arrays in shape [[x0,y0,z0], [x1,y1,z1],...] and surface is a vtk vtkDataSetSurfaceFilter object
    """

    streams = make_streamlines(vlsvfile, streamline_seeds, seeds_n, seeds_x0, seeds_range, dl, iterations)
    magnetopause = make_magnetopause(streams, end_x, x_point_n, sector_n, ignore)
    vertices, faces = make_surface(magnetopause)
    surface = make_vtk_surface(vertices, faces)

    return vertices, surface


   
