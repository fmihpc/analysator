"""Script and functions for creating sidecar files with SDF/region/boundary tags of plasma regions.

    variables used to find regions/boundary regions:
    rho, temperature, beta, beta_star, 



"""

import analysator as pt
import numpy as np
import vtk
#from .analysator.scripts import shue


R_E = 6371000

### Signed distance field functions: ###


def vtkDelaunay3d_SDF(all_points, coordinates):
    """Gives a signed distance to a convex hull surface created from given coordinates.

        :param all_points: points ([x, y, z] coordinates in m) for which a signed distance to surface will be calculated
        :param coordinates: coordinates (array of [x, y, z]:s in m) that are used to make a surface.
        :returns: array of signed distances, negative sign: inside the surface, positive sign: outside the surface
    
    """

    points =vtk.vtkPoints()#.NewInstance()
    for i in range(len(coordinates)):
        points.InsertNextPoint(coordinates[i,0],coordinates[i,1],coordinates[i,2])
        polydata = vtk.vtkPolyData()
        polydata.SetPoints(points)
        polydata.Modified()

    # Delaunay convex hull
    delaunay_3d = vtk.vtkDelaunay3D()
    delaunay_3d.SetInputData(polydata)
    delaunay_3d.Update()


    data = delaunay_3d.GetOutput()

    surface = vtk.vtkDataSetSurfaceFilter()
    surface.SetInputData(data)
    surface.Update()

    implicitPolyDataDistance = vtk.vtkImplicitPolyDataDistance()
    implicitPolyDataDistance.SetInput(surface.GetOutput())

    convexhull_sdf = np.zeros(len(all_points))

    # SDFs
    for i,coord in enumerate(all_points):
        convexhull_sdf[i] = implicitPolyDataDistance.EvaluateFunction(coord)
    
    return convexhull_sdf




def treshold_mask(data_array, value):
    """Chooses and flags cells where variable values match those given. If value is given as float, relative tolerance can be given

	:param data_array: array to mask, e.g. output of f.read_variable(name="proton/vg_rho", cellids=-1)
	:param variable: str, variable name
	:param value: value/values to use for masking; a float or int for exact match, a (value, relative tolerance) tuple, or [min value, max value] list pair where either can be None for less than or eq./more than or eq. value
	:returns: 0/1 mask in same order as cellids, 1: variable value in array inside treshold values, 0: outside
    """

    if data_array is None or value is None:  # either variable isn't usable or treshold has not been given
        return None

    mask = np.zeros((len(data_array)))
	
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
        print("Treshold mask didn't match any values in array")
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

    print(flag_name+" written to file")



def make_region_flags(variable_dict, condition_dict, flag_type="01", mask=None):
    """Makes region flags

        example:
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
                print(key, "ignored")

    else: # search only masked area
        for key in condition_dict.keys():
            try: # if key in variables
                region = treshold_mask(variable_dict[key][mask], condition_dict[key])
                if region is not None: variable_flags.append(region)
            except:
                print(key, "ignored")


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
        bowshock_rho_SDF = vtkDelaunay3d_SDF(query_points, f.get_cell_coordinates(cellids[bowshock_rho_flags!=0]))
        return bowshock_rho_SDF
    
    else: # bowshock from user-defined variable and values 
        bowshock_rho_flags = make_region_flags(variable_dict, own_condition_dict, flag_type="01")
        bowshock_rho_SDF = vtkDelaunay3d_SDF(query_points, f.get_cell_coordinates(cellids[bowshock_rho_flags!=0]))
        dict_sdf = vtkDelaunay3d_SDF(query_points, f.get_cell_coordinates(cellids[bowshock_rho_flags!=0]))

        return dict_sdf
    

def magnetopause_SDF(f, datafile, vtpoutfile, variable_dict, query_points, method="beta_star_with_connectivity", own_variable_dict=None): # TODO: remove vlsvReader AND datafile name need (streamline magnetopause to only use f?)
    """Finds the magnetopause by making a convex hull of either streamlines, beta*, or user-defined variables, and returns signed distances from surface to query_points.
        Note: constructing the magnetopause using solar wind streamlines is slow and needs more memory than e.g. the beta*-method
    
        :param query_points: xyz-coordinates of all points where the SDF will be calculated ([x1, y1, z1], [x2, y2, z2], ...)
        :kword method: str, specifies the method used to find magnetopause, options: "beta_star", "streamlines", "shue", "dict". If "dict" is used, own_variable_dict dictionary must be given. If "shue", run needs to be specified (TODO kwarg)
        :kword own_variable_dict: when using "dict"-method, dictionary with string variable names as keys and treshold pairs as values to pass to treshold_mask()-function
    """

    if method != "streamlines": cellids =f.read_variable("CellID")

    if method == "streamlines": 

        [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
        seeds_x0=150e6
        dl=5e5
        iters = int(((seeds_x0-xmin)/dl)+100)
        sector_n = 36*2
        vertices, vtkSurface = pt.calculations.find_magnetopause_sw_streamline_3d(datafile, seeds_n=200, seeds_x0=seeds_x0, seeds_range=[-15*6371000, 15*6371000], 
                                                                               dl=dl, iterations=iters, end_x=xmin+15*6371000, x_point_n=100, sector_n=sector_n) 

        # make the magnetopause surface from vertice points
        np.random.shuffle(vertices) # helps Delaunay triangulation
        Delaunay_SDF = vtkDelaunay3d_SDF(query_points, vertices)

        writer = vtk.vtkXMLPolyDataWriter() 
        writer.SetInputConnection(vtkSurface.GetOutputPort())
        writer.SetFileName(vtpoutfile)
        writer.Write()
        return Delaunay_SDF

    elif method == "beta_star":
        # magnetopause from beta_star only
        condition_dict = {"beta_star": [0.4, 0.5]} # FIC: [0.4, 0.5]) # EGE: [0.9, 1.0]) # max 0.6 in FHA to not take flyaways from outside magnetopause
        mpause_flags = make_region_flags(variable_dict, condition_dict, flag_type="01")
        contour_coords = f.get_cell_coordinates(cellids[mpause_flags!=0])
    
        # make a convex hull surface with vtk's Delaunay
        magnetopause_sdf = vtkDelaunay3d_SDF(query_points, contour_coords)
        return magnetopause_sdf
    

    elif method == "beta_star_with_connectivity":
        # magnetopause from beta_star, with connectivity if possible
        try:
            connectivity_region = treshold_mask(variable_dict["connection"], 0)
            betastar_region = treshold_mask(variable_dict["beta_star"], [0.6, 0.7])
            magnetosphere_proper =  np.where((connectivity_region==1) | (betastar_region==1), 1, 0)
            contour_coords = f.get_cell_coordinates(cellids[magnetosphere_proper==1])
        except:
            print("using field line connectivity for magnetosphere did not work, using only beta*")
            condition_dict = {"beta_star": [0.5, 0.6]} # FIC: [0.4, 0.5]) # EGE: [0.9, 1.0]) # max 0.6 in FHA to not take flyaways from outside magnetopause
            mpause_flags = make_region_flags(variable_dict, condition_dict, flag_type="01")
            contour_coords = f.get_cell_coordinates(cellids[mpause_flags!=0])
        
        # make a convex hull surface with vtk's Delaunay
        magnetopause_sdf = vtkDelaunay3d_SDF(query_points, contour_coords)
        return magnetopause_sdf
    
    elif method == "dict":
        # same method as beta* but from user-defined variables as initial contour
        flags = make_region_flags(variable_dict, own_variable_dict, flag_type="01")
        treshold_coords = f.get_cell_coordinates(cellids[flags!=0])
        dict_sdf = vtkDelaunay3d_SDF(query_points, treshold_coords)

        return dict_sdf
    
    elif method == "shue":
        return 0
        theta = np.linspace(0, 2*np.pi/3 , 200) # magnetotail length decided here by trial and error: [0, 5*np.pi/6] ~ -350e6 m, [0, 2*np.pi/3] ~ -100e6 m in EGE
        r, __, __ = shue.f_shue(theta, run='EGI') # EGI~EGE, for runs not in shue.py B_z, n_p, and v_sw need to be specified and run=None

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
        shue_SDF = vtkDelaunay3d_SDF(query_points, coords)

        return shue_SDF
    
    else:
        print("Magnetopause method not recognized. Use one of the options: \"beta_star\", \"beta_star_with_connectivity\", \"streamlines\", \"shue\", \"dict\"")

def RegionFlags(datafile, outfilen, vtpoufile, ignore_boundaries=True,
                magnetopause_method="beta_star", magnetopause_dict=None):
    
    """Creates a sidecar .vlsv file with flagged cells for regions and boundaries in near-Earth plasma environment. 
        Region flags (start with flag_, flags are fractions of filled conditions or 1/0): magnetosheath, magnetosphere, cusps, lobe_N, lobe_S, central_plasma_sheet, PSBL
        Boundary signed distance flags (start with SDF_, flags are signed distances to boundary in m with inside being negative distance): magnetopause, bowshock

    :param datafile: vlsv file name (and path)
    :param outdir: sidecar file save directory path
    :param outfilen: sidecar file name
    :kword ignore_boundaries: True: do not take cells in the inner/outer boundaries of the simulation into account when looking for regions 
    :kword magnetopause_method: default "beta_star", 
    """


    f = pt.vlsvfile.VlsvReader(file_name=datafile)
    cellids =f.read_variable("CellID")
    [xmin, ymin, zmin, xmax, ymax, zmax] = f.get_spatial_mesh_extent()
    all_points = f.get_cell_coordinates(cellids) #centre_points_of_cells(f)
    cellIDdict = f.get_cellid_locations()



    ## writer ##
    writer = pt.vlsvfile.VlsvWriter(f, outfilen)
    writer.copy_variables_list(f, ["CellID"])
    writer.copy_variables(f, varlist=["proton/vg_rho" , "vg_beta_star", "vg_temperature", "vg_b_vol", "vg_J", "vg_beta", "vg_connection", "vg_boundarytype"])


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


    def errormsg(varstr): print("{} could not be read, will be ignored".format(varstr))

    for varname, filevarname in varnames.items():
        if isinstance(filevarname, str): # no operator
            try: 
                variables[varname] = f.read_variable(name=filevarname, cellids=-1)
            except:
                errormsg(filevarname)
        else:
            try: 
                variables[varname] = f.read_variable(name=filevarname[0], cellids=-1, operator=filevarname[1])
            except:
                errormsg(filevarname)

    ## cellid testing
    #testmesh = np.where(cellids < 10000, 1, 0)
    #write_flags(writer, testmesh, 'testmesh')



    # upstream point
    upstream_point = [xmax-10*R_E,0.0,0.0]
    upstream_cellid = f.get_cellid(upstream_point)
    upstream_index = cellIDdict[upstream_cellid]

    ## MAGNETOPAUSE ##
    magnetopause = magnetopause_SDF(f, datafile, vtpoufile, variables, all_points, method=magnetopause_method, own_variable_dict=magnetopause_dict)
    write_flags(writer, magnetopause, 'SDF_magnetopause')
    write_flags(writer, np.where(np.abs(magnetopause) < 5e6, 1, 0), "flag_magnetopause")

    # save some magnetopause values for later
    magnetopause_density = np.mean(variables["density"][np.abs(magnetopause) < 5e6])
    print(f"{magnetopause_density=}")
    magnetopause_temperature = np.mean(variables["temperature"][np.abs(magnetopause) < 5e6])
    print(f"{magnetopause_temperature=}")
    
    ## MAGNETOSPHERE ##
    # magnetosphere from magentopause SDF
    magnetosphere = np.where(magnetopause<0, 1, 0)
    write_flags(writer, magnetosphere, 'flag_magnetosphere_convex')

    # magnetospshere from beta* and field line connectivity if possible # TODO
    try: 
        connectivity_region = treshold_mask(variables["connection"], 0)
        betastar_region = treshold_mask(variables["beta_star"], [0.01, 0.5])
        magnetosphere_proper =  np.where((connectivity_region==1) | (betastar_region==1), 1, 0)
        write_flags(writer, magnetosphere_proper, 'flag_magnetosphere')
    except:
        print("Non-convex magnetosphere could not be made")

    ## BOW SHOCK ##
    # bow shock from rho
    bowshock = bowshock_SDF(f, variables, all_points)
    write_flags(writer, bowshock, 'SDF_bowshock')
    write_flags(writer, np.where(np.abs(bowshock) < 5e6, 1, 0), "flag_bowshock")

    # magnetosphere+magnetosheath -area
    inside_bowshock = np.where(bowshock<0, 1, 0)

    ## MAGNETOSHEATH ##
    # magnetosheath from bow shock-magnetosphere difference
    magnetosheath_flags = np.where((inside_bowshock & 1-magnetosphere), 1, 0)
    write_flags(writer, magnetosheath_flags, 'flag_magnetosheath')

    # save magentosheath density and temperature for further use
    magnetosheath_density = np.mean(variables["density"][magnetosheath_flags == 1])
    print(f"{magnetosheath_density=}")
    magnetosheath_temperature = np.mean(variables["temperature"][magnetosheath_flags == 1])
    print(f"{magnetosheath_temperature=}")

    ## UPSTREAM ##
    # upstream from !bowshock
    write_flags(writer, 1-inside_bowshock, 'flag_upstream')

    write_flags(writer, inside_bowshock, 'flag_inside_bowshock')



    ## INNER MAGNETOSPHERE REGIONS ##

    if ignore_boundaries:
        noBoundaries = np.where(f.read_variable(name=boundaryname, cellids=-1) == 1, 1, 0) # boundarytype 1: not a boundary
        mask_inBowshock= np.where(((inside_bowshock == 1) & (noBoundaries == 1)), 1, 0).astype(bool)  # only search inner regions from inside the magnetosheath and magnetosphere
        mask_inMagnetosphere = np.where(((magnetosphere == 1) & (noBoundaries == 1)), 1, 0).astype(bool) # 

    else:
        mask_inBowshock = inside_bowshock.astype(bool)  # only search inner regions from inside the magnetosheath and magnetosphere
        mask_inMagnetosphere = magnetosphere.astype(bool) # 

    #print("upstream B:", variables["B_magnitude"][upstream_index])

    # cusps
    try:
        cusp_conditions = {"density": [variables["density"][upstream_index], None],
                        #"beta_star": [0.1, None],
                        "connection": [0.0, 2.5], # either closed, or open-closed/closed-open
                        "B_magnitude":[2*variables["B_magnitude"][upstream_index], None],
                        "J_magnitude": [variables["J_magnitude"][upstream_index], None]
                        }
    except:
        cusp_conditions = {"density": [variables["density"][upstream_index], None],
                        #"beta_star": [0.1, None],
                        "connection": [0.0, 2.5], # either closed, or open-closed/closed-open
                        "B_magnitude":[2*variables["B_magnitude"][upstream_index], None]
                        }

    cusp_flags = make_region_flags(variables, cusp_conditions, flag_type="fraction", mask=mask_inMagnetosphere)
    write_flags(writer, cusp_flags, 'flag_cusps', mask_inMagnetosphere)


    # magnetotail lobes 
    lobes_conditions = {"beta": [None, 0.1], # Koskinen: 0.003
                       "connection": [0.5, 2.5],
                       "density": [None, variables["density"][upstream_index]], # Koskinen: 1e-8
                        "temperature": [None, 3.5e6], # Koskinen: 3.5e6 K
                        }
    lobes_flags = make_region_flags(variables, lobes_conditions, flag_type="fraction", mask=mask_inBowshock)
    write_flags(writer, lobes_flags, 'flag_lobes', mask_inBowshock)

    # lobes slightly other way
    lobe_N_conditions = {"beta": [None, 0.1], 
                       "connection": [0.5, 2.5],
                       "B_x":[0, None],
                       "B_magnitude":[None, 10*variables["B_magnitude"][upstream_index]]
                       }
    lobe_N_flags = make_region_flags(variables, lobe_N_conditions, flag_type="fraction", mask=mask_inBowshock)
    write_flags(writer, lobe_N_flags, 'flag_lobe_N', mask_inBowshock)

    lobe_S_conditions = {"beta": [None, 0.1], 
                       "connection": [0.5, 2.5],
                       "B_x":[None, 0],
                       "B_magnitude":[None, 10*variables["B_magnitude"][upstream_index]]
                        }
    lobe_S_flags = make_region_flags(variables, lobe_S_conditions, flag_type="fraction", mask=mask_inBowshock)
    write_flags(writer, lobe_S_flags, 'flag_lobe_S', mask_inBowshock)

    # lobe density from median densities?
    lobes_mask = np.where((lobes_flags > 0.9), 1, 0).astype(bool)
    lobes_allcells = mask_inBowshock.astype(float)
    lobes_allcells[mask_inBowshock] = lobes_mask
    lobe_density = np.mean(variables["density"][lobes_allcells>0.9])
    print(f"{lobe_density=}")

    # Central plasma sheet
    central_plasma_sheet_conditions = {#"density": [None, variables["density"][upstream_index]],
                                        "density": [None, 1e6], # Wolf intro 
                                        #"density" : [1e7, None], # Koskinen: 3e-7
                                        #"density": [lobe_density, None],
                                       "beta": [1.0, None], # Koskinen: 6
                                       #"connection": 0, # only closed-closed
                                       #"temperature": [2*variables["temperature"][upstream_index], None]#,
                                       "temperature": [2e6, None], #  5e7 from Koskinen
                                       #"B_magnitude":[10*variables["B_magnitude"][upstream_index], None]
                                       #"J_magnitude": [2*variables["J_magnitude"][upstream_index], None],
                                       "J_magnitude": [1e-9, None],
                                       }
    
    central_plasma_sheet_flags = make_region_flags(variables, central_plasma_sheet_conditions,flag_type="fraction", mask=mask_inMagnetosphere)
    write_flags(writer, central_plasma_sheet_flags, 'flag_central_plasma_sheet', mask_inMagnetosphere)

    
    # Plasma sheet boundary layer (PSBL)
    #PSBL_conditions = {#"density": (1e7, 1.0), # Koskinen: 1e5
    #                    "density": [None, 1e7],
                        #"temperature": [0.5e6, 1e6],#(1e7, 2.0), # from Koskinen
                       #"temperature": [2e6, None],
    #                   "beta": [0.1, 1.0],
                       #"B_magnitude":[10*variables["B_magnitude"][upstream_index], None]
    #                   "J_magnitude": [2*variables["J_magnitude"][upstream_index], None],
    #                   }
    

    #PSBL_flags = make_region_flags(variables, PSBL_conditions,flag_type="fraction", mask=mask_inMagnetosphere)
    #write_flags(writer, PSBL_flags, "flag_PSBL", mask_inMagnetosphere)


    # Low-Latitude boundary layer (LLBL)
    #LLBL_conditions = {"density": [magnetopause_density, magnetosheath_density],
    #                   "temperature": [magnetosheath_temperature, magnetopause_temperature]
    #                   }
    #LLBL_flags = make_region_flags(variables, LLBL_conditions,flag_type="fraction", mask=mask_inBowshock)
    #write_flags(writer, LLBL_flags, "flag_LLBL", mask_inBowshock)





    # High-Latitude boundary layer (HLBL)
    HLBL_conditions = { "density": [magnetopause_density, magnetosheath_density],
                        "temperature": [magnetosheath_temperature, magnetopause_temperature],
                        "beta": [1.0, None], 
                        "connection": 0, 
                        "J_magnitude": [2*variables["J_magnitude"][upstream_index], None],
                        }
    
    HLBL_flags = make_region_flags(variables, HLBL_conditions,flag_type="fraction", mask=mask_inBowshock)
    write_flags(writer, HLBL_flags, "flag_HLBL", mask_inBowshock)



def main():

    datafile = "/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/bulk1.0001400.vlsv"
    outfilen = "FHA_regions_t0001400.vlsv"
    vtp_outfilen = "FHA_regions_t0001400.vtp"


    for infile, outfile, vtp_outfile in zip(datafile, outfilen, vtp_outfilen):
        RegionFlags(infile, outfile, vtp_outfile, magnetopause_method="streamlines")


if __name__ == "__main__":

    main()
