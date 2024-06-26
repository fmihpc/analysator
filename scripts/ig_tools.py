'''
    some tools for doing basic computations over the ionosphere mesh

    Note that the ionosphere is implemented as a triangular (refined Fibonacci) mesh, approximating a sphere

    Some terminology: "element" = triangular face of the mesh
                      "node" = corner of an element

    could be added as methods for VlsvReader objects in vlsvReader.py
'''

import numpy as np

def ionosphere_mesh_area(f):
    '''
       Inputs:
           f: VlsvReader object

       Output: 2D numpy array, areas of the triangular elements [m^2]
     
       this function could be added to analysator:master vlsvReader
    '''
    n = f.get_ionosphere_node_coords()       # nodes: shape (n_nodes, 3) vertices
    c = f.get_ionosphere_element_corners()   # corners of elements: indices integers 0-(n_nodes-1), shape (n_elements, 3)
    p = n[c,:]                               # shape(n_elements, 3, 3)   1st index is the element, 2nd index is triangle corners, 3rd index is x-y-z position
    r1 = p[:,1,:] - p[:,0,:]
    r2 = p[:,2,:] - p[:,0,:]
    areas = np.linalg.norm(np.cross(r1, r2), axis = 1) / 2.
    # checked: sum of triangle areas is near the expected area 4*pi*R^2 for a sphere:
    # ( np.sum(areas) - (np.pi * 4 ) * (R_EARTH + 100000.)**2 ) / np.sum(areas) ~ 1
    return areas


def get_ig_r(f):
    '''
       Inputs:
           f: VlsvReader object
       Output:
           barycenters of ionospheric elements. [n_elements, 3] numpy array

       Note: ionospheric mesh is at radius (R_EARTH + 100 km). See Ganse et al. (submitted 2024)

       this function could be added to analysator:master vlsvReader
    '''
    n = f.get_ionosphere_node_coords()          # Nodes. shape (n_nodes, 3)
    ec = f.get_ionosphere_element_corners()     # Element corners. shape (n_elements, 3)
    ig_r = np.zeros(np.array(ec).shape)
    for i in range(ig_r.shape[0]):
        ig_r[i,:] = (n[ec[i,0], :] + n[ec[i,1], :] + n[ec[i,2], :]) / 3  #barycenter, aka centroid
    return ig_r


def interp_ig_nodes_to_elements(f, var):
    '''
       Inputs:
           f: VlsvReader object
           var: Ionospheric variable defined on ionospheric grid corners
       Output:
        Barycentric interpolation, i.e. evaluated at the ig_ grid elements (triangles)
        
        this is here implemented as a linear interpolation to the barycenter,
        which turns out to be just adding 3 corner values and dividing by 3)

        assume var.shape[0] = n_nodes. Still works if var.ndim>1? (CHECK THIS)
    '''
    c = f.get_ionosphere_element_corners()                        # (n_elements, 3) node indices
    return np.nansum(var[c], axis = -1) / 3.


def _normalize(vec):
    '''
        (private) helper function, normalizes a multidimensinonal array of vectors
        assume [...., 3] array
    '''
    return vec / np.linalg.norm(vec, axis = -1)[:, np.newaxis]
    

def ig_E(f):
    '''
       Input:
           f: VlsvReader object
       Output:
           (in-plane) ionospheric electric field [V/m], derived from ionospheric potential 'ig_potential'

        The in-plane electric field has 2 degrees of freedom (orientation and magnitude),
        uniquely matching the potential at 2 corners of the triangular element.

        The algorithm:
            1. WLOG let one corner have potential V0=0. Let r1, r2 be position vectors for the other corners, at potentials V1 & V2.
            2. Construct an orthonormal basis {b1, b3} that spans the triangular element, where b1 || r1
            3. calculate the unique electric field E1 || b1, that gives the correct potential V1 at position r1
            4. infer the potential V1_B at the location r1_B (see figure)
            5. calculate the unique electric field E3 || b3, that is required to match the potentials at r1_B and r2.
            6. Calculate the total electric field E = E1 + E3
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
    ig_potential= f.read_variable('ig_potential')   # shape (n_nodes)
    c = f.get_ionosphere_element_corners()          # Element corners. shape (n_elements, 3)
    n = f.get_ionosphere_node_coords()       # Nodes. shape (n_nodes, 3)
    p = n[c,:]                               # shape(n_elements, 3, 3)   1st index is the element, 2nd identifies triangle corners, 3rd is the x-y-z position     
    r1 = p[:,1,:] - p[:,0,:]
    r2 = p[:,2,:] - p[:,0,:]
    r_shape = r1.shape
    d1 = np.repeat(np.linalg.norm(r1, axis = 1), 3).reshape(r_shape)   # r1 distance (repeated 3 times for convenience)
    # define some orthogonal basis vectors b1 and b3 (where b1 || r1) that span the triangular element:
    b1 = _normalize(r1)
    r1_B = b1 * np.repeat(np.nansum(r2 * b1, axis = 1), 3).reshape(r_shape)
    r3 = r2 - r1_B
    b3 = _normalize(r3)
    d3 = np.repeat(np.linalg.norm(r3, axis = 1), 3).reshape(r_shape)
    # electric field E = -grad V:
    ig_potential_c = ig_potential[c]
    V1 = np.repeat(ig_potential_c[:,1] - ig_potential_c[:,0], 3).reshape(r_shape)  # potential diff. btw. nodes 1 and 0 of a given triangular face
    V2 = np.repeat(ig_potential_c[:,2] - ig_potential_c[:,0], 3).reshape(r_shape) 
    E1 = -b1 * V1 / d1         # (E1 || b1)
    V1_B =  np.repeat(np.nansum(-E1 * r1_B, axis = 1), 3).reshape(r_shape)
    E3 = -b3 * (V2 - V1_B) / d3         # (E_1 || b1)
    E = E1 + E3
    # checked: V1 === np.nansum(-E * r1, axis = 1),  and   V2 === np.nansum(-E * r2, axis = 1)
    return E



def ig_inplanecurrent_from_ig_E(f):
    '''
        Input:
           f: VlsvReader object
        Output:
            calculate height-integrated in-plane current [A/m]for ig_inplanecurrent from the ionospheric electric field ig_E

        This is probably only needed to reconstruct the ig_inplanecurrent for .vlsv files
        that happen to be missing this variable (as in t=501-1000 in run FHA)
    '''
    E = ig_E(f)
    ig_sigmah = np.repeat(interp_ig_nodes_to_elements(f, f.read_variable('ig_sigmah')), 3).reshape(E.shape) 
    ig_sigmap = np.repeat(interp_ig_nodes_to_elements(f, f.read_variable('ig_sigmap')), 3).reshape(E.shape) 
    ig_r = get_ig_r(f)
    ig_r_hat = _normalize(ig_r)
    # v SIGN v   ( approximate b_hat = r_hat in southern hemisphere, and b_hat = -r_hat in northern hemisphere):
    ig_b_hat = -ig_r_hat * np.repeat(np.sign(ig_r[:,2]), 3).reshape(ig_r_hat.shape) 
    return ig_sigmap * E - ig_sigmah * np.cross(E, ig_b_hat)  # see 'Currents-in-the-ionosphere.pdf, eq. 5.7'



def ig_inplanecurrent(f):
    '''
       Wrapper function. Returns height-integrated in-plane current,
       whether or not the 'ig_inplanecurrent' variable is saved in the file.

       This is probably only needed to reconstruct the ig_inplanecurrent for .vlsv files
       that happen to be missing this variable (as in t=501-1000 in run FHA)

       Input:
           f: VlsvReader object
       Output: height-integrated in-plane current [A/m],
    '''
    try:
        ipc = f.read_variable('ig_inplanecurrent')
    except ValueError:  # if 'ig_inplanecurrent' not in vlsvReader object
        ipc = ig_inplanecurrent_from_ig_E(f)
    return ipc
