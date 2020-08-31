import pytools as pt
import plot_colormap #gives error but doesn't work without 
import sys, os, socket
import numpy as np
from operator import attrgetter

import yt
import matplotlib.pylab as plt
from yt.visualization.api import Streamlines

from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap


'''
running:
-takes two (or zero) commad line arguments: how many inner streamlines to ignore and from how many next to count the magnetopause position

Some issues: 
-Nose part unstable and slows the whole thing down, also requires lots of streamlines to work
-Bulk file change is hard and must be done by hand
-Everything else must also be done by changing the code by hand
'''



class Point:
    '''
    a 2D point
    '''

    def __init__(self, x, z):
        self.x = x
        self.z = z

        #radius and angle in regard to origo
        self.r, self.phi = np.array(cartesian_to_polar(x,z))

        #radius to x-axis
        self.r_x = np.sqrt(z**2)

        #radius to z-axis
        self.r_z = np.sqrt(x**2)



class Slice: #slice of a plane

    def __init__(self, plane, constant, pos):
        self.plane = plane # 'YZ' or 'XY'
        self.constant = constant # x or z -coordinate depending on the plane
        self.pos = pos # above or below the plane, above = 0, below = 1
        self.points = [] #contains points in slice


    def add_point(self, p):
        self.points.append(p)


    def get_mean_point(self, ignore, count):
    #get mean of points in slice s ignoring first 'ignore' points
    #and counting the next 'count' slices
    
        #if self.plane == 'XY':
        #    ignore = 1
        #    count = 2

        if len(self.points) < (ignore+count): #if there are not enough points in slice just ignore it
            return 0

        
        #sort points in slice by their radius
        if self.plane == 'YZ': 
            pnts = sorted(self.points, key=attrgetter('r_x'))
        elif self.plane == 'XY': 
            pnts = sorted(self.points, key=attrgetter('r_z')) 


        m = [] #list to save the radius of points that count
        i = 0 #index of point under consideration
        for p in pnts:
            if i < ignore: #too close
                i = i + 1
                continue
            elif i < (ignore+count): #points that count
                if self.plane == 'YZ':
                    m.append(p.r_x)
                elif self.plane == 'XY':
                    m.append(p.r_z)
                
            else: #went over
                break
            i = i + 1

        #get mean of the collected radii
        mean_r = np.mean(m)

        #below or above
        if self.pos == 0: r = mean_r
        if self.pos == 1: r = -mean_r

        #make (x,z)-coordinates of the mean point
        if self.plane == 'YZ': 
            coord_mean = (self.constant, r)
        if self.plane == 'XY':
            coord_mean = (r, self.constant)
        return coord_mean


class Plane:

    def __init__(self, plane, constant):
        self.plane = plane # 'XY' or 'YZ'
        self.constant = constant # z or x -coordiante depending on the plane
        self.slices = []

        if self.plane == 'YZ':
            self.slices.append(Slice(self.plane, self.constant, 0)) #above
            self.slices.append(Slice(self.plane, self.constant, 1)) #below
        elif self.plane == 'XY':
            self.slices.append(Slice(self.plane, self.constant, 0)) #just above
    

    def add_point(self, p): # add a point to the right slice in plane
        
        if self.plane == 'YZ':
            if p.phi > 0 and p.phi <= 180:
                ind = 0
            else:
                ind = 1 
        elif self.plane == 'XY':
            ind = 0
    
        sl = self.slices[ind]
        sl.add_point(p)




def get_bulk(): #fetch bulk file name

    #BCQ: (bulks from 0001500 to 0002500)

    run = 'BCQ'
    num = 1550
    num = str(num)

    fileLocation="/wrk/group/spacephysics/vlasiator/2D/"+run+"/bulk/"
    fileN = "bulk.000"+num+".vlsv"

    f = fileLocation+fileN
    return f, run



def to_Re(m): #meters to Re
    Re = 6371000
    return m/Re

def cartesian_to_polar(x, z): #cartesian coordinates to polar
    r = np.sqrt(x**2 + z**2)
    phi = np.degrees(np.arctan2(z, x))

    
    return(r, phi)

def polar_to_cartesian(r, phi): #polar coordinates to cartesian
    phi = np.radians(phi)

    z = r * np.sin(phi)
    x = r * np.cos(phi)

    return(x, z)


def make_streamlines():
    boxre=[0,0]
                            
    # bulk file
    name = get_bulk()[0]
    f = pt.vlsvfile.VlsvReader(file_name=name)

    #get box coordinates from data
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    [xsize, ysize, zsize] = f.get_spatial_mesh_size()
    simext=[xmin,xmax,ymin,ymax,zmin,zmax]
    sizes=np.array([xsize,ysize,zsize])

    simext = list(map(to_Re, simext))

    #set box coordnates
    if len(boxre)==4:
        boxcoords=boxre
    else:
        boxcoords=list(simext)

   
    # If box extents were provided manually, truncate to simulation extents
    boxcoords[0] = max(boxcoords[0],simext[0]) #xmax
    boxcoords[1] = min(boxcoords[1],simext[1]) #xmin
    boxcoords[2] = max(boxcoords[2],simext[2]) #ymax
    boxcoords[3] = min(boxcoords[3],simext[3]) #ymin
    boxcoords[4] = max(boxcoords[4],simext[4]) #zmax
    boxcoords[5] = min(boxcoords[5],simext[5]) #zmin
    
    cellids = f.read_variable("CellID")
    
        
    #Read the data from vlsv-file
    Vx = f.read_variable("v", operator="x")
    Vz = f.read_variable("v", operator="z")

    #Re-shape variable data
    Vxs=Vx[np.argsort(cellids)].reshape(f.get_spatial_mesh_size(), order="F")
    Vys = np.zeros_like(Vxs)
    Vzs=Vz[np.argsort(cellids)].reshape(f.get_spatial_mesh_size(), order="F")
    

    data=dict(Vx=Vxs,Vy=Vys,Vz=Vzs)
    
    
    #Create streamline seeds (starting points for streamlines)
    streamline_seeds = []
        
    #range: np.arange(from, to, step)
    for i in np.arange(-10.0, 10.0, 0.01): # 2000 seeds is nice, 400 is on the low end
        streamline_seeds.append([20, 0, i])

            
    streamline_seeds = np.array(streamline_seeds)


    #dataset in yt-form
    yt_dataset = yt.frontends.stream.load_uniform_grid(
        data,
        sizes,
        bbox=np.array([[boxcoords[0], boxcoords[1]],
                       [boxcoords[2],boxcoords[3]],
                       [boxcoords[4],boxcoords[5]]]),
        periodicity=[True, True, True]) # Has to be forced...                                 

    #data, seeds, dictionary positions, lenght of lines
    streamlines_pos = yt.visualization.api.Streamlines(yt_dataset, streamline_seeds,
                                                       "Vx", "Vy", "Vz", length=60, direction=1)

    #where the magic happens
    streamlines_pos.integrate_through_volume()
        

    return streamlines_pos #returns yt streamlines object


def get_poynting_flux(coordinate_points):

    #input to meters
    Re = 6371000
    coordinate_points = [[coord*Re for coord in arr] for arr in coordinate_points]

    #data file
    name = get_bulk()[0]
    f = pt.vlsvfile.VlsvReader(file_name=name)

    #get middle points of next-to-next coordinates
    middle_points = []
    normals = []
    for i,j in zip(np.arange(0, len(coordinate_points)-1), np.arange(1, len(coordinate_points))):
        A = coordinate_points[i]
        B = coordinate_points[j]

        middle = np.mean([A, B], axis=0)
        middle_points.append(middle)

        # normal vectors of lines between next-to-next coordinate points:
        n = [-(B[2]-A[2]), 0, (B[0]-A[0])] #the outwards-facing surface vector
        normal_vector = n/np.linalg.norm(n) #normalization
        normals.append(normal_vector)

    

    middle = [[coord/Re for coord in arr] for arr in middle_points]


    #get interpolated poynting vectors at magnetopause coordinate middle points
    S = f.read_interpolated_variable("poynting", middle_points)

    #Pounting flux (S.n)
    Poynting_flux = []
    Poynting_flux_normal =[]
    for s, n in zip(S, normals):
        P_flux = np.dot(s,n)
        Poynting_flux.append(P_flux)
        Poynting_flux_normal.append(P_flux*n)

   

    
    return Poynting_flux, Poynting_flux_normal, middle



def interpolate(YTstream, points, axis):
    '''
    :kword YTstream:    A single streamline in yt Streamline -form
    :kword points:      Array of some axis points where coordinates of streamline passing a plane are wanted
    :kword axis:        x or z depending on what axis the points are on
    '''

    streamline = YTstream.positions.to_value()
    arr = np.array(streamline)
    
    #maybe needs to be sorted???
    xp = arr[:,0] #x-coordinates of streamline
    #yp = arr[:,1] #y-coordinates
    zp = arr[:,2] #z-coordinates

    #reverse arrays so interpolation works
    xp = xp[::-1]
    zp = zp[::-1]

    valid = True
    #interpolate missing coordinate at points
    if axis == 'x':
        x_points = points
        z_points = np.interp(points, xp, zp, left=float('inf'), right=float('inf'))
        for z in z_points:
            if z == float('inf'): #if z-coordinate is inf at some point scrap the whole streamline from data
                valid = False

        
        #gather interpolated points to new stream coordinates
        stream_points = [x_points, z_points]

        return stream_points, valid

    elif axis == 'z': 
        x_points = []
        z_points = []

        for z0 in points:
            x_coords = []
            for i in range(0, len(xp)-1):
                if (zp[i] < z0 and zp[i+1] > z0) or (zp[i] > z0 and zp[i+1] < z0):
                    x_coords.append(xp[i] + (z0-zp[i])*(xp[i+1]-xp[i])/(zp[i+1]-zp[i]))

            if len(x_coords) > 0:
                clean_x_coords = [x for x in x_coords if x > 5]
                if len(clean_x_coords) > 0:
                    x_points.append(min(clean_x_coords))
                else:
                    x_points.append(float('inf'))
            else:
                x_points.append(float('inf'))

        stream_points = [x_points, points]
        return stream_points, True

    else:
        print('axis must be x or z!')
        exit() 


def get_magnetopause(ignore, count):
    
    streams = make_streamlines()


    #the tail of the magnetopause

    #wanted points in x-axis (YZ-plane places)
    x_points = np.array(np.arange(-20, 7, 0.2))

    #planes
    YZplanes = []
    for x in x_points:
        YZplanes.append(Plane('YZ', x))
    
    #go through all streamlines one by one
    for s in np.arange(0, len(streams.streamlines)):
        
        YTstream = streams.path(s)

        #interpolate z-values of x_points            
        stream_points, valid = interpolate(YTstream, x_points, 'x')
        
        if valid:
            #add interpolated points to yz-planes
            for i in range(0, len(x_points)): #for every YZ-plane

                if stream_points[0][i] != x_points[i]: #check that values for x match
                    print("Something's very wrong")
                    exit()

                p = Point(stream_points[0][i], stream_points[1][i])
                YZplanes[i].add_point(p)
        

    #list to save the magnetopause body coordinates by slice (above/below x-axis)
    pause_coords = [[],[]]

    #for every yz-plane
    for p in YZplanes:

        #for every slice get mean point
        for i in range(0, len(p.slices)):
            s = p.slices[i] #slice

            mean_point = s.get_mean_point(ignore, count)
            
            if mean_point != 0:
                pause_coords[i].append(mean_point)
    
    #get the ends of the body lines for the nose
    above_end = pause_coords[0][-1]
    below_end = pause_coords[i][-1]




    #the nose

    #XY-plane places
    dz=0.2
    z_points = np.array(np.arange(below_end[1]+dz, above_end[1]-dz, dz))

    XYplanes = []
    for z in z_points:
        XYplanes.append(Plane('XY', z))

    
    #go through all streamlines one by one again
    for s in np.arange(0, len(streams.streamlines)):
        YTstream = streams.path(s)

        #linear interpolation of the points in the nose
        tip_points, valid = interpolate(YTstream, z_points, 'z')
        
        if valid:
            #add interpolated points to XY-planes
            for i in range(0, len(z_points)):
                if tip_points[0][i] != float('inf'): 
                    p = Point(tip_points[0][i], tip_points[1][i])
                    XYplanes[i].add_point(p)

    #construct the nose
    
    nose = []
    
    for p in XYplanes: #for every plane
        for s in p.slices: #for every slice

            mean_point = s.get_mean_point(ignore, count)
            if mean_point != 0:
                nose.append(mean_point)

    nose = np.array(nose)


    
    #make a continous line of the magnetopause
    #structure: above, nose, below
    coords = []
    for i in 0,1: 
        sl = pause_coords[i] #body coordinates

        if len(sl) == 0: #if slice is empty (shouldn't be)
            continue

        if i == 1: #make the nose in second loop
            nose = nose[::-1]
            for x, z in zip(nose[:,0], nose[:,1]):
                coords.append([x, 0, z])

            #also reverse the below-line so it starts from where the nose ended
            sl = sl[::-1] 
        
        
        sl = np.array(sl)
        for x, z in zip(sl[:,0], sl[:,1]):
            coords.append([x, 0, z])

    
    
    return (streams, coords)  #returns streamlines, array of magnetopause positions by sloce and by point, poynting vector ends in 2d-array 


def plot_all(ax, XmeshXY=None, YmeshXY=None, pass_maps=None): #external function for plot_colormap

    #what to plot, boolean
    streamlines = False #streamlines
    magnetopause = True #magnetopause line
    poynting_flux_vectors = False # S.n * normal vectors
    poynting_flux = False #poynting flux as colours


    #command line arguments (how many inner lines to ignore and count)
    if len(sys.argv) < 2:
        ignore, count = 3, 5
    else: 
        ignore = int(sys.argv[1])
        count = int(sys.argv[2])


    streams, coords = get_magnetopause(ignore, count)
    

    if streamlines:
        for s in np.arange(0, len(streams.streamlines)):
            stream_pos = streams.streamlines[s]
            ax.plot(stream_pos[:,0], stream_pos[:,2], color='#C4C4C4', lw=0.1)



    if magnetopause:
        coords = np.array(coords)
        ax.plot(coords[:,0], coords[:,2], 'k', lw = 0.6)



    if poynting_flux_vectors:
        poynting, vectors, middle_coords = get_poynting_flux(coords)

        coords = np.array(middle_coords)
        p = np.array(vectors)
        p = p*10

        ax.quiver(coords[:,0], coords[:,2], p[:,0], p[:,2], width=0.001)

        

    if poynting_flux:
        poynting = get_poynting_flux(coords)[0]

        coords = np.array(coords)
        p = np.array(poynting)*100

        maxval = abs(max(p, key=abs))

        x = coords[:,0]
        z = coords[:,2]
       
        points = np.array([x, z]).T.reshape(-1,1,2)
        segments = np.concatenate([points[:-1],points[1:]], axis=1)

        lc = LineCollection(segments, cmap='PiYG', linewidth=1.5, norm=plt.Normalize(vmax=maxval, vmin=-maxval))
        lc.set_array(p)
        plt.gca().add_collection(lc)



def main():

    inputfile, run_name = get_bulk()
    
    plot_colormap.plot_colormap(
        filename=inputfile,
        #step=j,
        nooverwrite=None,
        var="rho",
        boxre = [-25,20,-20,20],
        title=None,
        draw=None, usesci=True, Earth=None,
        external = plot_all,
        run=run_name,
        #vectors='B', vectordensity=500, vectorsize=0.5
        #streamlines="B", streamlinedensity=1.5, streamlinethick=0.7,
        )  
    

if __name__ == "__main__":
    main()




