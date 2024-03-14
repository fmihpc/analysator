# some tools for dealing with the ionosphere mesh
# could be added as methods for VlsvReader objects in vlsvReader.py

import numpy as np


# duplicate of ionosphere_mesh_area() as it appears in biot_savart.py
def ionosphere_mesh_area(f):
    # this function could be added to analysator:master vlsvReader
    n = f.get_ionosphere_node_coords()       # nodes: shape (21568, 3) vertices
    c = f.get_ionosphere_element_corners()   # corners of elements: indices integers 0-21567, shape (43132, 3)
    p = n[c,:]                               # shape(43132, 3, 3)   first index is the element, second index identifies corners of the triangle, third index is the x-y-z position
    r1 = p[:,1,:] - p[:,0,:]
    r2 = p[:,2,:] - p[:,0,:]
    areas = np.linalg.norm(np.cross(r1, r2), axis = 1) / 2.     # use cross product to calculate areas
    # checked: sum of triangle areas is near the expected area 4*pi*R^2 for a sphere
    # ( np.sum(areas) - (np.pi * 4 ) * (R_EARTH + 100000.)**2 ) / np.sum(areas) 
    return areas


# duplicate of get_ig_r() as it appears in biot_savart.py
def get_ig_r(f):
    # Calculate barycenters of ionospheric elements. Could add this to master:analysator vlsvReader
    # ionospheric mesh is at radius (R_EARTH + 100 km), see Urs's ionosphere writeup
    n = f.get_ionosphere_node_coords()          # node = vertex of the triangular mesh
    ec = f.get_ionosphere_element_corners()     # (Element Corners), where element = triangular face
    #ig_r = np.zeros(ec.shape)
    ig_r = np.zeros(np.array(ec).shape)
    for i in range(ig_r.shape[0]):
        ig_r[i,:] = (n[ec[i,0], :] + n[ec[i,1], :] + n[ec[i,2], :]) / 3  #barycenter, aka centroid
    return ig_r    # [n_elements, 3]


'''
def read_variable_at_ig_nodes(f, ig_var):
    #    read an ionospheric variable, evaluated at the ig_ grid nodes (triangle corners)
    #    this may require upsampling if the variable is only given at the elements (triangles)
    return f.get_ionosphere_node_coords() * 0. # (n_nodes, 3) DUMMY
'''

def interp_ig_nodes_to_elements(reader, var):
    '''
        interpolate, evaluated at the ig_ grid elements (triangles)
        this is here implemented as a linear interpolation to the barycenter
        (which turns out to be just adding 3 corner values and dividing by 3)

        assume var.shape[0] = n_nodes. Still works if var.ndim>1? (CHECK THIS)
    '''
    c = reader.get_ionosphere_element_corners()                        # (n_elements, 3) node indices
    return np.nansum(var[c], axis = -1) / 3.


def normalize(vec):
    '''
        assume [...., 3] array
    '''
    return vec / np.repeat(np.linalg.norm(vec, axis = -1), 3).reshape(vec.shape)


#add this as a data reducer
def ig_E(reader):
    '''
        calculate (in-plane) ionospheric electric field from ig_potential
               .
              ^^\
     r2,(V2) / | 
            /  | 
           /   | r3 (d3), b3, E3
          /    |
         /    _|     \
        .____|_._____>.  r1 (d1), b1, E1, (V1)
    V0=0    r1_B, V1_B

       E = E1 + E3
    '''
    ig_potential= reader.read_variable('ig_potential')   # shape (n_nodes)
    c = reader.get_ionosphere_element_corners()                        # (n_elements, 3)
    n = reader.get_ionosphere_node_coords()       # nodes: shape (21568, 3) vertices
    p = n[c,:]                               # shape(43132, 3, 3)   1st index is the element, 2nd identifies triangle corners, 3rd is the x-y-z position     
    r1 = p[:,1,:] - p[:,0,:]
    r2 = p[:,2,:] - p[:,0,:]
    r_shape = r1.shape
    d1 = np.repeat(np.linalg.norm(r1, axis = 1), 3).reshape(r_shape)   # r1 distance (repeated 3 times for convenience)
    # define some orthogonal basis vectors b1 and b3 (where b1 || r1) that span the triangular element:
    b1 = normalize(r1)
    r1_B = b1 * np.repeat(np.nansum(r2 * b1, axis = 1), 3).reshape(r_shape)
    r3 = r2 - r1_B
    b3 = normalize(r3)
    d3 = np.repeat(np.linalg.norm(r3, axis = 1), 3).reshape(r_shape)
    # electric field E = -grad V:
    ig_potential_c = ig_potential[c]
    V1 = np.repeat(ig_potential_c[:,1] - ig_potential_c[:,0], 3).reshape(r_shape)  # potential diff. btw. nodes 1 and 0 of a given triangular face
    V2 = np.repeat(ig_potential_c[:,2] - ig_potential_c[:,0], 3).reshape(r_shape) 
    E1 = -b1 * V1 / d1         # (E_1 || b1)
    V1_B =  np.repeat(np.nansum(-E1 * r1_B, axis = 1), 3).reshape(r_shape)
    E3 = -b3 * (V2 - V1_B) / d3         # (E_1 || b1)
    E = E1 + E3
    # checked: V1 === np.nansum(-E * r1, axis = 1),  and   V2 === np.nansum(-E * r2, axis = 1)
    return E



def ig_inplanecurrent_from_ig_E(reader):
    '''
        (re-)calculate for ig_inplanecurrent from the ionospheric electric field ig_E
        This is probably only needed to reconstruct the ig_inplanecurrent for .vlsv files
        that happen to be missing this variable (as in t=501-1000 in run FHA)
    '''
    E = ig_E(reader)
    ig_sigmah = np.repeat(interp_ig_nodes_to_elements(reader, reader.read_variable('ig_sigmah')), 3).reshape(E.shape) 
    ig_sigmap = np.repeat(interp_ig_nodes_to_elements(reader, reader.read_variable('ig_sigmap')), 3).reshape(E.shape) 
    ig_r = get_ig_r(reader)
    ig_r_hat = normalize(ig_r)
    # v SIGN v   ( approximate b_hat = r_hat in southern hemisphere, and b_hat = -r_hat in northern hemisphere):
    ig_b_hat = -ig_r_hat * np.repeat(np.sign(ig_r[:,2]), 3).reshape(ig_r_hat.shape) 
    return ig_sigmap * E - ig_sigmah * np.cross(E, ig_b_hat)  # see 'Currents-in-the-ionosphere.pdf, eq. 5.7'



def ig_inplanecurrent(reader):
    try:
        ipc = reader.read_variable('ig_inplanecurrent')
    except ValueError:  # if 'ig_inplanecurrent' not in vlsvReader object
        ipc = ig_inplanecurrent_from_ig_E(reader)
    return ipc
