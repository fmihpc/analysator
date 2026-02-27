"""Functions for finding the magnetopause from .vlsv data

    Example use with own arguments passed to streamline function:
    
    .. code-block:: python

        datafile = "bulkfile.vlsv"
        vtpoutfilen = "magnetopause.vtp"
        vlsvoutfilen = "magnetopause.vlsv"

        SW_args = {"seeds_n":300, "seeds_x0":150e6, "seeds_range":[-5*6371000, 5*6371000], 
                "dl":5e5, "iterations":int(1500), "end_x":-50*6371000, "x_point_n":100, "sector_n":36*2} 

        surface, SDF = magnetopause(datafile,
                                method="streamlines",
                                return_SDF=True,
                                return_surface=True,
                                method_args=SW_args) 

        write_vtk_surface_to_file(surface, vtpoutfilen)
        write_SDF_to_file(SDF, datafile, vlsvoutfilen)

        
    Example use for 0.0-0.4 beta* magnetopause with tetrahedrons with sides longer than 1Re excluded (aplha = 1Re):

    .. code-block:: python

        datafile = "bulkfile.vlsv"
        vtpoutfilen = "magnetopause.vtp"

        surface, __ = magnetopause(datafile,
                                method="beta_star_with_connectivity", 
                                beta_star_range=[0.0, 0.4],
                                Delaunay_alpha=1*6371000,
                                return_SDF=False,
                                return_surface=True) 

        write_vtk_surface_to_file(surface, vtpoutfilen)


"""

import analysator as pt
import numpy as np
import vtk
from scripts import shue, regions
import logging


R_E = 6371000
    

def write_vtk_surface_to_file(vtkSurface, outfilen):
    writer = vtk.vtkXMLPolyDataWriter() 
    writer.SetInputConnection(vtkSurface.GetOutputPort())
    writer.SetFileName(outfilen)
    writer.Write()
    logging.info("wrote ", outfilen)

def write_SDF_to_file(SDF, datafilen, outfilen):
    f = pt.vlsvfile.VlsvReader(file_name=datafilen)
    writer = pt.vlsvfile.VlsvWriter(f, outfilen)
    writer.copy_variables_list(f, ["CellID"])
    writer.write_variable_info(pt.calculations.VariableInfo(SDF, "SDF_magnetopause", "-", latex="",latexunits=""),"SpatialGrid",1)
    logging.info("wrote ", outfilen)


      

def magnetopause(datafilen, method="beta_star_with_connectivity", own_tresholds=None, return_surface=True, return_SDF=True, SDF_points=None, Delaunay_alpha=None, beta_star_range=[0.4, 0.5], method_args={}): # TODO: separate streamline suface and vtkDelaunay3d surface in streamline method
    """Finds the magnetopause using the specified method. Surface is constructed using vtk's Delaunay3d triangulation which results in a convex hull if no Delaunay_alpha is given.
        Returns vtk.vtkDataSetSurfaceFilter object and/or signed distances (negative -> inside magnetopause) (=SDF) to all cells
        Note that using alpha for Delaunay might make SDF different from expected inside the magnetosphere, especially if surface is constructed with points not everywhere in the magnetosphere (e.g. beta* 0.4-0.5) or if simulation grid size is larger than alpha

    :param datafilen: a .vlsv bulk file name (and path)
    :kword method: str, default "beta_star_with_connectivity", other options "beta_star", "streamlines", "shue", "dict"
    :kword own_tresholds: if using method "dict", a dictionary with conditions for the magnetopause/magnetosphere must be given where key is data name in datafile and value is treshold used (see treshold_mask())
    :kword return_surface: True/False, return vtkDataSetSurfaceFilter object
    :kword return_SDF: True/False, return array of distances in m to SDF_points in point input order, negative distance inside the surface
    :kword SDF_points: optionally give array of own points to calculate signed distances to. If not given, distances will be to cell centres in the order of f.read_variable("CellID") output
    :kword Delaunay_alpha: alpha (float) to give to vtkDelaunay3d, None -> convex hull, alpha=__: surface egdes longer than __ will be excluded
    :kword beta_star_range: [min, max] treshold rage to use with methods "beta_star" and "beta_star_with_connectivity"
    :kword method_args: dict of keyword arguments to be passed down to external functions (for streamlines and shue)
    :returns: vtkDataSetSurfaceFilter object of convex hull or alpha shape if return_surface=True, signed distance field of convex hull or alpha shape of magnetopause if return_SDF=True
    """

    f = pt.vlsvfile.VlsvReader(file_name=datafilen)
    cellids = f.read_variable("CellID")
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    query_points = f.get_cell_coordinates(cellids) # for SDF, all centre points of cells

    if method == "streamlines": 

        if "streamline_seeds" in method_args or "seeds_n" in method_args:
            vertices, manual_vtkSurface = pt.calculations.find_magnetopause_sw_streamline_3d(datafilen, **method_args)
        else:
            seeds_x0=150e6
            dl=5e5
            iters = int(((seeds_x0-xmin)/dl)+100)
            sector_n = 36*3
            vertices, manual_vtkSurface = pt.calculations.find_magnetopause_sw_streamline_3d(datafilen, seeds_n=300, seeds_x0=seeds_x0, seeds_range=[-5*6371000, 5*6371000], 
                                                                                    dl=dl, iterations=iters, end_x=xmin+10*6371000, x_point_n=200, sector_n=sector_n) 

        # make the magnetopause surface from vertice points
        np.random.shuffle(vertices) # helps Delaunay triangulation
        vtkSurface, SDF = regions.vtkDelaunay3d_SDF(query_points, vertices, Delaunay_alpha)

    elif method == "beta_star":
        # magnetopause from beta_star only
        mpause_flags = regions.treshold_mask(f.read_variable("vg_beta_star"), beta_star_range)
        contour_coords = f.get_cell_coordinates(cellids[mpause_flags!=0])
    
        # make a convex hull surface with vtk's Delaunay
        np.random.shuffle(contour_coords)
        vtkSurface, SDF = regions.vtkDelaunay3d_SDF(query_points, contour_coords, Delaunay_alpha)

    

    elif method == "beta_star_with_connectivity":
        # magnetopause from beta_star, with connectivity if possible
        betastar_region = regions.treshold_mask(f.read_variable("vg_beta_star"), beta_star_range)
        try:
            connectivity_region = regions.treshold_mask(f.read_variable("vg_connection"), 0) # closed-closed magnetic field lines
            magnetosphere_proper =  np.where((connectivity_region==1) | (betastar_region==1), 1, 0)
            contour_coords = f.get_cell_coordinates(cellids[magnetosphere_proper==1])
            np.save("pointcloud.npy", contour_coords)
        except:
            logging.warning("using field line connectivity for magnetosphere did not work, using only beta*")
            #condition_dict = {"beta_star": [0.5, 0.6]} # FIC: [0.4, 0.5]) # EGE: [0.9, 1.0]) # max 0.6 in FHA to not take flyaways from outside magnetopause
            mpause_flags = np.where(betastar_region==1, 1, 0)
            contour_coords = f.get_cell_coordinates(cellids[mpause_flags!=0])
        
        # make a convex hull surface with vtk's Delaunay
        vtkSurface, SDF = regions.vtkDelaunay3d_SDF(query_points, contour_coords, Delaunay_alpha)

    #elif method == "beta_star_with_fieldlines": # either incredibly slow or does not work, don't use without fixing #TODO proprer measure of actual field line backwall point
    #    # magnetopause from beta_star, with field lines connecting to back wall if possible
    #    betastar_region = regions.treshold_mask(f.read_variable("vg_beta_star"), beta_star_range)
    #    try:
    #        fw_region = regions.treshold_mask(f.read_variable("vg_connection_coordinates_fw")[:,0], [None, xmin+5e7])
    #        bw_region = regions.treshold_mask(f.read_variable("vg_connection_coordinates_bw")[:,0], [None, xmin+5e7])
    #        #connectivity_region = regions.treshold_mask(f.read_variable("vg_connection"), 0) # closed-closed magnetic field lines
    #        magnetosphere_proper =  np.where(((fw_region==1) | (betastar_region==1) | (bw_region==1)), 1, 0)
    #        contour_coords = f.get_cell_coordinates(cellids[magnetosphere_proper==1])
    #    except:
    #        logging.warning("using field lines for magnetosphere did not work, using only beta*")
    #        #condition_dict = {"beta_star": [0.5, 0.6]} # FIC: [0.4, 0.5]) # EGE: [0.9, 1.0]) # max 0.6 in FHA to not take flyaways from outside magnetopause
    #        mpause_flags = np.where(betastar_region==1, 1, 0)
    #        contour_coords = f.get_cell_coordinates(cellids[mpause_flags!=0])
    #    
    #    # make a convex hull surface with vtk's Delaunay
    #    vtkSurface, SDF = regions.vtkDelaunay3d_SDF(query_points, contour_coords, Delaunay_alpha)


    
    elif method == "dict":
        variable_dict = {}
        for key in own_tresholds:
            variable_dict[key] = f.read_variable(key, cellids=-1)

        # same method as beta* but from user-defined variables as initial contour
        flags = regions.make_region_flags(variable_dict, own_tresholds, flag_type="01")
        treshold_coords = f.get_cell_coordinates(cellids[flags!=0])
        vtkSurface, SDF = regions.vtkDelaunay3d_SDF(query_points, treshold_coords, Delaunay_alpha)

    
    elif method == "shue":
        # note: might not be correct but should produce something, should recheck the projection coordinates

        theta = np.linspace(0, 2*np.pi/3 , 200) # magnetotail length decided here by trial and error: [0, 5*np.pi/6] ~ -350e6 m, [0, 2*np.pi/3] ~ -100e6 m in EGE
        if "run" in method_args or "B_z" in method_args:
             r, __, __ = shue.f_shue(theta, **method_args)
        else:
            r, __, __ = shue.f_shue(theta, B_z = -10, n_p = 1, v_sw = 750) # for runs not in shue.py B_z, n_p, and v_sw need to be specified and run=None

        # 2d one-sided magnetopause
        xs = r*np.cos(theta)
        ys = r*np.sin(theta)

        # 3d projection for complete magnetopause
        psis = np.linspace(0, 2*np.pi, 100)
        coords = np.zeros((len(theta)*len(psis), 3))
        i=0
        for x,y in zip(xs,ys): 
            for psi in psis:
                coords[i] = np.array([x, y*np.sin(psi), y*np.cos(psi)])
                i += 1

        coords = coords*R_E

        # surface and SDF
        np.random.shuffle(coords) # helps Delaunay triangulation
        vtkSurface, SDF = regions.vtkDelaunay3d_SDF(query_points, coords) # alpha does nothing here


    
    else:
        logging.error("Magnetopause method not recognized. Use one of the options: \"beta_star\", \"beta_star_with_connectivity\", \"streamlines\", \"shue\", \"dict\"")
        exit()

    if return_surface and return_SDF:
        return vtkSurface, SDF
    elif return_surface:
        return vtkSurface, None
    elif return_SDF:
        return None, SDF
    
    # beta* + B connectivity:
    #connectivity_region = regions.treshold_mask(f.read_variable("vg_connection"), 0)
    #betastar_region = regions.treshold_mask(f.read_variable("beta_star"), [0.0, 0.5])
    #magnetosphere_proper =  np.where((connectivity_region==1) | (betastar_region==1), 1, 0)


def main():

    datafile = "/wrk-vakka/group/spacephysics/vlasiator/3D/EGE/bulk/bulk.0002000.vlsv"
    vtpoutfilen = "EGE_magnetopause_t2000.vtp"
    vlsvoutfilen = "EGE_magnetopause_t2000.vlsv"

    surface, SDF = magnetopause(datafile,
                            method="beta_star", 
                            beta_star_range=[0.9, 1.0],
                            return_SDF=True,
                            return_surface=True) 

    write_vtk_surface_to_file(surface, vtpoutfilen)
    write_SDF_to_file(SDF, datafile, vlsvoutfilen)


if __name__ == "__main__":

    main()

