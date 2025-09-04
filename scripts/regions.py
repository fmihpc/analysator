"""Script and functions for creating sidecar files with SDF/region/boundary tags of plasma regions.

    Usage example where bow shock conditions are given to replace default conditions:

    .. code-block:: python
    
        datafile = "bulkfile.vlsv"
        outfilen = "regions.vlsv"
        RegionFlags(datafile, outfilen, regions=["all"],
                    region_conditions = {"bowshock": {"density": [2e6, None]}}

"""

import analysator as pt
import numpy as np
import vtk
from scripts import magnetopause
import logging

R_E = 6371000

def vtkDelaunay3d_SDF(query_points, coordinates, alpha=None):
    """Gives a signed distance to a convex hull or alpha shape surface created from given coordinates.
        Note: if using alpha, SDF might not work, especially if the point cloud used for Delaunay does not fully cover e.g. lobes

        :param all_points: points ([x, y, z] coordinates in m) for which a signed distance to surface will be calculated
        :param coordinates: coordinates (array of [x, y, z]:s in m) that are used to make a surface.
        :kword alpha: alpha to be given to vtkDelaunay3D (e.g. R_E), removes surface edges that have length more than alpha so that the resulting surface is not convex. None -> convex hull
        :returns: vtkDataSetSurfaceFilter() object, array of signed distances (negative sign: inside the surface, positive sign: outside the surface)
    
    """

    points = vtk.vtkPoints()#.NewInstance()
    for i in range(len(coordinates)):
        points.InsertNextPoint(coordinates[i,0],coordinates[i,1],coordinates[i,2])
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.Modified()

    # Delaunay convex hull
    delaunay_3d = vtk.vtkDelaunay3D()
    delaunay_3d.SetInputData(polydata)
    if alpha is not None:
        delaunay_3d.SetAlpha(alpha)
        delaunay_3d.SetAlphaTets(1)
        delaunay_3d.SetAlphaLines(0)
        delaunay_3d.SetAlphaTris(0)
        delaunay_3d.SetAlphaVerts(0)
    delaunay_3d.Update()


    data = delaunay_3d.GetOutput()

    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(data)
    surface.Update()

    implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
    implicitPolyDataDistance.SetInput(surface.GetOutput())

    convexhull_sdf = np.zeros(len(query_points))

    # SDFs
    for i,coord in enumerate(query_points):
        convexhull_sdf[i] = implicitPolyDataDistance.EvaluateFunction(coord)
    
    return surface, convexhull_sdf


def vtkSDF(query_points, dualgrid):
    ''' Obtain the SDF in relation to a waterproof volumetric grid "dualgrid"
    '''

    # print(dualgrid)
    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(dualgrid)
    surface.Update()

    # print(surface.GetOutput())


    implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
    implicitPolyDataDistance.SetInput(surface.GetOutput())
    import sys
    convexhull_sdf = np.zeros(len(query_points))
    for i,coord in enumerate(query_points):
        convexhull_sdf[i] = implicitPolyDataDistance.EvaluateFunction(coord)
        # sys.exit()
    
    return surface.GetOutputDataObject(0), convexhull_sdf


def treshold_mask(data_array, value):
    """Chooses and flags cells where variable values match those given. If value is given as float, relative tolerance can be given

	:param data_array: array to mask, e.g. output of f.read_variable(name="proton/vg_rho", cellids=-1)
	:param variable: str, variable name
	:param value: value/values to use for masking; a float or int for exact match, a (value, relative tolerance) tuple, or [min value, max value] list pair where either can be None for less than or eq./more than or eq. value
	:returns: 0/1 mask in same order as data_array, 1: variable value in array inside treshold values, 0: outside
    """

    if (data_array is None) or (value is None) or (np.isnan(data_array[0])):  # either variable isn't usable or treshold has not been given
        return None

    #mask = np.zeros((len(data_array)))
	
    if isinstance(value, float) or isinstance(value, int): # single value, exact
        mask = np.where(np.isclose(data_array, value), 1, 0)
    elif isinstance(value, tuple): # single value with tolerance
        mask = np.where(np.isclose(data_array, value[0], rtol=value[1]), 1, 0)
    else:
        if value[0] is None:
            mask = np.where((data_array <= value[1]), 1, 0) # anywhere where value is less than or equal to
        elif value[1] is None:
            mask = np.where((value[0] <= data_array), 1, 0) # anywhere where value is more than or equal to
        else:
            mask = np.where((value[0] <= data_array) & (data_array <= value[1]), 1, 0) # min/max

    if np.sum(mask[mask>0]) == 0:
        logging.warning("Treshold mask didn't match any values in array")
        return None

    return mask


def box_mask(f, coordpoints=None, cells=None, marginal=[150e6, 50e6, 50e6, 50e6, 50e6, 50e6]):
    """Crops simulation box for calculations, output flags outside of cropped box will be 0

        :param f: a VlsvReader
        :kword coordpoints: coordinate points of cells to be masked
        :kword cells: cellIDs to be masked
        :kword marginal: 6-length list of wanted marginal lengths from mesh edges in meters [negx, negy, negz, posx, posy, posz]
        :returns: Boolean mask in input order of coordinates/cells, -1 if no coorpoints or cells were given and mask could not be done
    """
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    
    def is_inside(coords):
        res = np.zeros((len(coords)), dtype=bool)

        for i,coord in enumerate(coords):
            if np.any((coord[0] < xmin+marginal[0], coord[0] > xmax-marginal[3], 
                      coord[1] < ymin+marginal[1], coord[1] > ymax-marginal[4], 
                      coord[2] < zmin+marginal[2], coord[2] > zmax-marginal[5])):
                res[i] = False
            else:
                res[i] =True
        return res
    
    if coordpoints is not None:
        return is_inside(coordpoints)
    
    elif cells is not None:
        coords = f.get_cell_coordinates(cells)
        return is_inside(coords)

    return -1

def write_flags(writer, flags, flag_name, mask=None):
    """Writes flags into .vlsv-file with writer

        :param writer: a vlsvwriter
        :param flags: an array of flags in order of cellids
        :param flag_name: string, name of the flag variable
        :kword mask: if flags are made from cropped cells, give mask used for cropping
    """

    if mask is not None: # flags exist only in cropped cells
        new_flags = mask.astype(float) 
        new_flags[mask] = flags
        writer.write_variable_info(pt.calculations.VariableInfo(new_flags, flag_name, "-", latex="",latexunits=""),"SpatialGrid",1)
    
    else:
        #writer.write(flags,flag_name,'VARIABLE','SpatialGrid')
        writer.write_variable_info(pt.calculations.VariableInfo(flags, flag_name, "-", latex="",latexunits=""),"SpatialGrid",1)

    logging.info(flag_name+" written to file")



def make_region_flags(variable_dict, condition_dict, flag_type="01", mask=None):
    """Makes region flags

        example:
        
        .. code-block:: python
        
            make_region_flags(variable_dict = {"density": f.read_variable(name="proton/vg_rho", cellids=-1)}, condition_dict = {"density": [None, 1e7]})

        results in flag array in shape of read_variable output where cells where "proton/vg_rho" is less or equal to 1e7 are 1 and others are 0
        
    
        :param variable_dict: dictionary containing names and data arrays of variables, names must match those in condition_dict
        :param condition_dict: dictionary containing names and conditions for variable data in variable_dict
        :kword flag_type: default "01", optionally "fraction". "01" gives flags 1: all variables fill all conditions, 0: everything else; "fraction" flags are rounded fractions of varibles that fulfill conditions (e.g. flag 0.5 means one of two conditions was filled)
        :kword mask: deafault None, optionally boolean array to use for variable dict arrays to search for a region from subset of cells (cells False in mask will be automatically flagged as 0)
    
    """

    variable_flags = []

    if mask is None: # search all cells
        for key in condition_dict.keys():

            try: # if key in variables
                region = treshold_mask(variable_dict[key], condition_dict[key])
                if region is not None: variable_flags.append(region)
            except:
                logging.warning(key, "ignored")

    else: # search only masked area
        for key in condition_dict.keys():
            try: # if key in variables
                region = treshold_mask(variable_dict[key][mask], condition_dict[key])
                if region is not None: variable_flags.append(region)
            except:
                logging.warning(key, "ignored")


    variable_flags = np.array(variable_flags)

    if len(variable_flags)==0: # no cells that fullfill conditions, hope that density is a key
        if mask is None:
            return np.zeros((len(variable_dict["density"])))
        else:
            return np.zeros((len(variable_dict["density"][mask])))
        

    if flag_type == "01": # 1 only if all viable variables are 1, 0 otherwise
        flags = np.where(np.sum(variable_flags, axis=0)==len(variable_flags), 1, 0)
        return flags
    
    elif flag_type == "fraction": # rounded fraction of viable variables that are 1
        flags = np.sum(variable_flags, axis=0)/len(variable_flags)
        return flags.round(decimals=2)
    


def bowshock_SDF(f, variable_dict, query_points, own_condition_dict=None):
    """Finds the bow shock by making a convex hull of either default variables/values based on upstream values or a user-defined variable/value dictionary. Returns signed distances from the convex hull to given query_points.
    
        :param f: a VlsvReader
        :param variable_dict: dictionary containing names and variable arrays needed
        :param query_points: xyz-coordinates of all points where the SDF will be calculated ([x1, y1, z1], [x2, y2, z2], ...)
        :kword own_condition_dict: optional, dictionary with string variable names as keys and tresholds as values to pass to treshold_mask()-function
    """

    cellids =f.read_variable("CellID")

    if own_condition_dict is None:
    # bow shock from upstream rho
        # upstream point
        upstream_point = [150e6, 0.0, 0.0]
        upstream_cellid = f.get_cellid(upstream_point)
        upstream_rho = f.read_variable(name="proton/vg_rho", cellids=upstream_cellid)

        bowshock_conditions = {"density": (1.5*upstream_rho, 0.1)}
        bowshock_rho_flags = make_region_flags(variable_dict, bowshock_conditions, flag_type="01")
        __, bowshock_rho_SDF = vtkDelaunay3d_SDF(query_points, f.get_cell_coordinates(cellids[bowshock_rho_flags!=0]))
        return bowshock_rho_SDF
    
    else: # bowshock from user-defined variable and values 
        bowshock_dict_flags = make_region_flags(variable_dict, own_condition_dict, flag_type="01")
        __, dict_sdf = vtkDelaunay3d_SDF(query_points, f.get_cell_coordinates(cellids[bowshock_dict_flags!=0]))

        return dict_sdf
    


def RegionFlags(datafile, outfilen, regions=["all"], ignore_boundaries=True, region_flag_type="01", magnetopause_kwargs={}, region_conditions={}):
    """Creates a sidecar .vlsv file with flagged cells for regions and boundaries in near-Earth plasma environment. 
        Region flags (named flag_region, flags are fractions of filled conditions or 1/0): magnetosheath, magnetosphere, cusps, lobe_N, lobe_S, central_plasma_sheet
        Boundary signed distance flags (named "SDF_boundary", flags are signed distances to boundary in m with inside being negative distance): magnetopause, bowshock

        possible regions: "all", "boundaries" (magnetopause, bow shock), "large_areas" (boundaries + upstream, magnetosheath, magnetosphere), "magnetosphere", "bowshock",
            "cusps", "lobes", "central_plasma_sheet"

        Note that different runs may need different tresholds for region parameters and region accuracy should be verified visually
        Beta* convex hull is most likely the best magnetopause method for regions, for just magnetopause with different options use magnetopause.py

        :param datafile: .vlsv bulk file name (and path)
        :param outfilen: sidecar .vlsv file name (and path)
        :kword ignore_boundaries: True: do not take cells in the inner/outer boundaries of the simulation into account when looking for regions (applicable for cusps and CPS for now)
        :kword region_flag_type: "01" or "fraction", whether flags are binary (all conditions must be satisfied) or fractions of how many given conditions are met
        :kword magnetopause_method: default "beta_star", other options: see magnetopause.py, use alpha=None
        :kword region_conditions: optional dict where keys are str region names and values are condition dictionaries, for setting own conditions for bow shock, cusps, lobes or central plasma sheet
    """

    if "all" in regions:
        regions.extend(["magnetopause", "bowshock", "upstream", "magnetosheath", "magnetosphere", "cusps", "lobes", "central_plasma_sheet"])
    else:
        if "large_areas" in regions:
            regions.extend(["magnetopause", "bowshock", "upstream", "magnetosheath", "magnetosphere"])
        if "boundaries" in regions:
            regions.extend(["magnetopause", "bowshock"])
        if "cusps" in regions or "central_plasma_sheet" in regions:  # for cusps and central plasma sheet we need the magnetoshpere
            regions.extend(["magnetosphere"])
        if "magnetosphere" in regions or "magnetosheath" in regions:
            regions.extend(["magnetopause"])
        if "bowshock" in regions or "magnetosheath" in regions or "upstream" in regions:
            regions.extend(["bowshock"])



    f = pt.vlsvfile.VlsvReader(file_name=datafile)
    cellids =f.read_variable("CellID")
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    all_points = f.get_cell_coordinates(cellids)
    cellIDdict = f.get_cellid_locations()

    ## writer ##
    writer = pt.vlsvfile.VlsvWriter(f, outfilen)
    writer.copy_variables_list(f, ["CellID"])
    #writer.copy_variables(f, varlist=["proton/vg_rho" , "vg_beta_star", "vg_temperature", "vg_b_vol", "vg_J", "vg_beta", "vg_connection", "vg_boundarytype"])


    # variable data dictionary
    variables = {}
    
    if f.check_variable("proton/vg_rho"): # vg-grid

        # dict {varible call name : variable vlsvreader name or varible call name : (variable vlsvreader name, operator)} 
        varnames = {"density": "proton/vg_rho",
                    "temperature": "vg_temperature",
                    "beta": "vg_beta",
                    "beta_star": "vg_beta_star",
                    "B_magnitude": ("vg_B_vol", "magnitude"),
                    "B_x": ("vg_B_vol", "x"),
                    "B_y": ("vg_B_vol", "y"),
                    "B_z": ("vg_B_vol", "z"),
                    "connection": "vg_connection",
                    "J_magnitude": ("vg_J", "magnitude")}
        
        boundaryname = "vg_boundarytype"


    def errormsg(varstr): logging.warning("{} could not be read, will be ignored".format(varstr))

    for varname, filevarname in varnames.items():
        if isinstance(filevarname, str): # no operator
            try: 
                variables[varname] = f.read_variable(name=filevarname, cellids=-1)
            except:
                variables[varname] = np.full((len(cellids)), np.nan)
                errormsg(filevarname)
        else:
            try: 
                variables[varname] = f.read_variable(name=filevarname[0], cellids=-1, operator=filevarname[1])
            except:
                variables[varname] = np.full((len(cellids)), np.nan)
                errormsg(filevarname)



    #connectivity_region = treshold_mask(variables["connection"], 0)
    #betastar_region = treshold_mask(variables["beta_star"], [0.0, 0.5])
    #magnetosphere_proper =  np.where((connectivity_region==1) | (betastar_region==1), 1, 0)
    #write_flags(writer, magnetosphere_proper, 'flag_magnetosphere')
    #return 0


    # upstream point # this could probably be done better than just a random point in sw
    upstream_point = [xmax-10*R_E,0.0,0.0]
    upstream_cellid = f.get_cellid(upstream_point)
    upstream_index = cellIDdict[upstream_cellid]

    ## MAGNETOPAUSE ##
    if "magnetopause" in regions:
        if magnetopause_kwargs:
            __, magnetopause_SDF = magnetopause.magnetopause(datafile, **magnetopause_kwargs)
        else:
            __, magnetopause_SDF = magnetopause.magnetopause(datafile, method="beta_star_with_connectivity", Delaunay_alpha=2e6, )  # default magnetopause: beta*+ B connectivity convex hull 
        write_flags(writer, magnetopause_SDF, 'SDF_magnetopause')
        write_flags(writer, np.where(np.abs(magnetopause_SDF) < 5e6, 1, 0), "flag_magnetopause")

        # save some magnetopause values for later
        magnetopause_density = np.mean(variables["density"][np.abs(magnetopause_SDF) < 5e6])
        #print(f"{magnetopause_density=}")
        magnetopause_temperature = np.mean(variables["temperature"][np.abs(magnetopause_SDF) < 5e6])
        #print(f"{magnetopause_temperature=}")
    
    ## MAGNETOSPHERE ##
    # magnetosphere from magnetopause SDF
    if "magnetosphere" in regions:
        magnetosphere = np.where(magnetopause_SDF<0, 1, 0)
        write_flags(writer, magnetosphere, 'flag_magnetosphere')


    ## BOW SHOCK ## #TODO: similar kwargs system as magnetopause?
    # bow shock from rho 
    if "bowshock" in regions:
        if "bowshock" in region_conditions:
            bowshock = bowshock_SDF(f, variables, all_points, own_condition_dict=region_conditions["bowshock"])
        else:
            bowshock = bowshock_SDF(f, variables, all_points) # default upstream rho method, might fail with foreshock
        write_flags(writer, bowshock, 'SDF_bowshock')
        write_flags(writer, np.where(np.abs(bowshock) < 5e6, 1, 0), "flag_bowshock")

        # magnetosphere+magnetosheath -area
        inside_bowshock = np.where(bowshock<0, 1, 0)

    ## MAGNETOSHEATH ##
    if "magnetosheath" in regions:
        # magnetosheath from bow shock-magnetosphere difference
        magnetosheath_flags = np.where((inside_bowshock & 1-magnetosphere), 1, 0)
        write_flags(writer, magnetosheath_flags, 'flag_magnetosheath')

        # save magnetosheath density and temperature for further use
        #magnetosheath_density = np.mean(variables["density"][magnetosheath_flags == 1])
        #print(f"{magnetosheath_density=}")
        #magnetosheath_temperature = np.mean(variables["temperature"][magnetosheath_flags == 1])
        #print(f"{magnetosheath_temperature=}")

    ## UPSTREAM ##
    if "upstream" in regions:
        # upstream from !bowshock
        write_flags(writer, 1-inside_bowshock, 'flag_upstream')
        #write_flags(writer, inside_bowshock, 'flag_inside_bowshock')



    ## INNER MAGNETOSPHERE REGIONS ##

    if "magnetosphere" in regions:
        if ignore_boundaries:
            noBoundaries = np.where(f.read_variable(name=boundaryname, cellids=-1) == 1, 1, 0) # boundarytype 1: not a boundary
            #mask_inBowshock= np.where(((inside_bowshock == 1) & (noBoundaries == 1)), 1, 0).astype(bool)  # only search inner regions from inside the magnetosheath and magnetosphere
            mask_inMagnetosphere = np.where(((magnetosphere == 1) & (noBoundaries == 1)), 1, 0).astype(bool) # 
        else:
            #mask_inBowshock = inside_bowshock.astype(bool)  # only search inner regions from inside the magnetosheath and magnetosphere
            mask_inMagnetosphere = magnetosphere.astype(bool) # 

    #print("upstream B:", variables["B_magnitude"][upstream_index])

    # cusps
    if "cusps" in regions:
        if "cusps" in region_conditions:
            cusp_conditions = region_conditions["cusps"]
        else:
            cusp_conditions = {"density": [variables["density"][upstream_index], None],
                            #"beta_star": [0.1, None],
                            "connection": [0.0, 2.5], # either closed, or open-closed/closed-open
                            "B_magnitude":[2*variables["B_magnitude"][upstream_index], None],
                            "J_magnitude": [variables["J_magnitude"][upstream_index], None]
                            }

        cusp_flags = make_region_flags(variables, cusp_conditions, flag_type=region_flag_type, mask=mask_inMagnetosphere)
        write_flags(writer, cusp_flags, 'flag_cusps', mask_inMagnetosphere)

 
    # magnetotail lobes 
    if "lobes" in regions:
        if "lobes" in region_conditions:
            lobes_conditions = region_conditions["lobes"]
        else:
            lobes_conditions = {"beta": [None, 0.1], # Koskinen: 0.003
                                "connection": [0.5, 2.5],
                                "density": [None, variables["density"][upstream_index]], # Koskinen: 1e-8
                                "temperature": [None, 3.5e6], # Koskinen: 3.5e6 K
                                }

            # lobes slightly other way
            lobe_N_conditions = {"beta": [None, 0.1], 
                                "connection": [0.5, 2.5],
                                "B_x":[0, None],
                                "B_magnitude":[None, 10*variables["B_magnitude"][upstream_index]]
                                }
            lobe_N_flags = make_region_flags(variables, lobe_N_conditions, flag_type=region_flag_type)
            

            lobe_S_conditions = {"beta": [None, 0.1], 
                                "connection": [0.5, 2.5],
                                "B_x":[None, 0],
                                "B_magnitude":[None, 10*variables["B_magnitude"][upstream_index]]
                                }
            lobe_S_flags = make_region_flags(variables, lobe_S_conditions, flag_type=region_flag_type)

            write_flags(writer, lobe_N_flags, 'flag_lobe_N')
            write_flags(writer, lobe_S_flags, 'flag_lobe_S')

        lobes_flags = make_region_flags(variables, lobes_conditions, flag_type=region_flag_type)
        write_flags(writer, lobes_flags, 'flag_lobes')


        # lobe density from median densities?
        #lobes_mask = np.where((lobes_flags > 0.9), 1, 0).astype(bool)
        #lobe_density = np.mean(variables["density"][lobes_mask])#[lobes_allcells>0.9])
        #print(f"{lobe_density=}")

    # Central plasma sheet
    if "central_plasma_sheet" in regions:
        if "central_plasma_sheet" in region_conditions:
            central_plasma_sheet_conditions = region_conditions["central_plasma_sheet"]
        else:
            central_plasma_sheet_conditions = {"density": [None, 1e6], # Wolf intro ?should not be but seems to work
                                            "beta": [1.0, None], # Koskinen: 6
                                            "temperature": [2e6, None], #  5e7 from Koskinen
                                            "J_magnitude": [1e-9, None],
                                             }
        
        central_plasma_sheet_flags = make_region_flags(variables, central_plasma_sheet_conditions,flag_type=region_flag_type, mask=mask_inMagnetosphere)
        write_flags(writer, central_plasma_sheet_flags, 'flag_central_plasma_sheet', mask_inMagnetosphere)


    ## Other boundary layers, PSBL sometimes works
    
    # Plasma sheet boundary layer (PSBL)
    # if "PSBL" in regions or "all" in regions:
    #PSBL_conditions = {#"density": (1e7, 1.0), # Koskinen: 1e5
    #                    "density": [None, 1e7],
                        #"temperature": [0.5e6, 1e6],#(1e7, 2.0), # from Koskinen
                       #"temperature": [2e6, None],
    #                   "beta": [0.1, 1.0],
                       #"B_magnitude":[10*variables["B_magnitude"][upstream_index], None]
    #                   "J_magnitude": [2*variables["J_magnitude"][upstream_index], None],
    #                   }
    
def main():

    fileid = 1000
    datafile = "/wrk-vakka/group/spacephysics/vlasiator/3D/FID/bulk1/bulk1.{:07d}.vlsv".format(fileid) 
    outfilen = "/wrk-vakka/group/spacephysics/vlasiator/3D/FID/postprocessing/prototyping/FID_mpause_{:07d}.vlsv".format(fileid) 

    RegionFlags(datafile, outfilen, regions=["magnetopause"])


if __name__ == "__main__":

    main()
