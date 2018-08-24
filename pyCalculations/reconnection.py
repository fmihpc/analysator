import pytools as pt
import numpy as np
from scipy.signal import convolve2d
import matplotlib.pyplot as plt
from shapely import geometry
from numpy import linalg as LA


# FUNCTION FOR FINDING THE INTERSECTION POINTS OF TWO CONTOURS
def findIntersection(v1,v2):
  # modified from https://stackoverflow.com/questions/42838190/finding-intersection-of-two-contour-plots-in-python
  #p1 = contour1.collections[0].get_paths()[0]
  #v1 = p1.vertices

  #p2 = contour2.collections[0].get_paths()[0]
  #v2 = p2.vertices

  poly1 = geometry.LineString(v1)
  poly2 = geometry.LineString(v2)

  intersection = poly1.intersection(poly2)

  return intersection



def find_X_O(fluxfile, bulkfile, step):

    vlsvfile = pt.vlsvfile.VlsvReader(bulkfile);

    ## Open bulkfile, determine sizes ##
    x_cells=int(vlsvfile.get_spatial_mesh_size()[0])
    z_cells=int(vlsvfile.get_spatial_mesh_size()[2])
    xsize = vlsvfile.read_parameter("xcells_ini")
    xmax =  vlsvfile.read_parameter("xmax")
    xmin =  vlsvfile.read_parameter("xmin")
    zmin =  vlsvfile.read_parameter("zmin")
    zmax =  vlsvfile.read_parameter("zmax")
    dx = (xmax-xmin)/xsize 

    ## DEFINE ARRAYS FOR AXIS
    x_array=np.array(range(int(xmin), int(xmax), int(dx)))
    z_array=np.array(range(int(zmin), int(zmax), int(dx)))

    # Open input fluxfile
    flux_function = np.fromfile(fluxfile,dtype='double').reshape(z_cells,x_cells)
    flux_offset = float(step)*-0.3535*5e-9*(-7.5e5)

    # Smooth fluxfunction
    kernel_size=5
    fkernel = np.ones((kernel_size,kernel_size))/(kernel_size**2)
    flux_function= convolve2d(flux_function, fkernel, 'same')

    # read rho for plotting later
    B = vlsvfile.read_variable( name="B", operator="z")
    cellids= vlsvfile.read_variable('CellID')
    B = B[cellids.argsort()].reshape(z_cells,x_cells)

    # calculate gradient of flux function
    dfdx,dfdz=np.gradient(flux_function)

    ##calculate the 0 contours of df/dx and df/dz
    #pl.figure(1)
    contour1=plt.contour(x_array,z_array, dfdx, [0], colors=('b'))
    contour1_paths=contour1.collections[0].get_paths()
    contour2=plt.contour(x_array,z_array, dfdz, [0], colors=('g'))
    contour2_paths=contour2.collections[0].get_paths()

    x_coords=[]
    z_coords=[]
    # find the intersection points of the 
    for path1 in contour1_paths:
       for path2 in contour2_paths:
          if path1.intersects_path(path2)  and len(path1)>1 and len(path2)>1:
             intersection=findIntersection(path1.vertices,path2.vertices)
             intersection_points=np.asarray(intersection)
             if len(intersection_points)>0:
                
                if len(intersection_points.shape)==1:
                   x_coords.append(intersection_points[0])
                   z_coords.append(intersection_points[1])
                else:
                   for i in range(len(intersection_points[:,0])):
                      x_coords.append(intersection_points[i,0])
                      z_coords.append(intersection_points[i,1])


    # DEFINE the type of the gradient(flux)=0 ##
    x_point_location=[]
    o_point_location=[]
    #minima_location=[]
    flux_function=flux_function.T

    for k in range(len(x_coords)):
       #cellid = 1+i+j*x_cells
       coords=[x_coords[k],0,z_coords[k]]
       cellid=vlsvfile.get_cellid(coords)
       i=int((cellid-1)%x_cells)
       j=(int(cellid)-1)/x_cells
       difference=[]
       
       ## the limist for i and j have the value 100, to save time, I have not been interested about x's and o's within 100 cell from boundaries.
       ## If you are then change the limits
       if i > 100 and j > 100 and i < x_cells-100 and j < z_cells-100:
          # Hessian matrix using central difference formulas for the second partial derivatives
          deltaPsi_xx= (flux_function[i+1,j]-2*flux_function[i,j]+flux_function[i-1,j])/dx**2
          deltaPsi_zz= (flux_function[i,j+1]-2*flux_function[i,j]+flux_function[i,j-1])/dx**2
          deltaPsi_xz= (flux_function[i+1,j+1]-flux_function[i+1,j-1]-flux_function[i-1,j+1]+flux_function[i-1,j-1])/(4*dx**2)
          
          Hessian = [[deltaPsi_xx, deltaPsi_xz], [deltaPsi_xz, deltaPsi_zz]]
          DetHess = deltaPsi_xx*deltaPsi_zz-deltaPsi_xz*deltaPsi_xz
          eigvals, eigvectors = LA.eig(Hessian)
          #coords.append(DetHess)
          #if sign_changes == 4:
          #if eigvals[0]*eigvals[1] <0 :
             
          if DetHess < 0:
             x_point_location.append(coords)           
                            
          ## NOTE if you want the o-points to be local maxima use deltaPsi_xx < 0, if you want them to be local minima use  deltaPsi_xx > 0                  
          #if DetHess > 0 and deltaPsi_xx > 0:
             #minima_location.append(coords)
          if DetHess > 0 and deltaPsi_xx < 0:
             o_point_location.append(coords)

    return x_point_location, o_point_location


