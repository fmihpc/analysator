# VLSV reader in Julia
#
# Hongyang Zhou, hyzhou@umich.edu 12/01/2020

include("vlsvvariables.jl")

using LightXML

"Mesh size information."
struct MeshInfo
   vxblocks::Int64
   vyblocks::Int64
   vzblocks::Int64
   vxblock_size::Int64
   vyblock_size::Int64
   vzblock_size::Int64
   vxmin::Float64
   vymin::Float64
   vzmin::Float64
   vxmax::Float64
   vymax::Float64
   vzmax::Float64
   dvx::Float64
   dvy::Float64
   dvz::Float64
end

"""
Variable metadata from the vlsv footer, including the unit of the 
variable as a regular string, the unit of the variable as a LaTeX-formatted 
string, the description of the variable as a LaTeX-formatted string, and the 
conversion factor to SI units as a string.
"""
struct VarInfo
   unit::String
   unitLaTeX::LaTeXString
   variableLaTeX::LaTeXString
   unitConversion::String
end

"Meta data declaration."
struct Meta
   name::AbstractString
   fid::IOStream
   footer::XMLElement
   fileindex_for_cellid::Vector{UInt64}
   cellIndex::Vector{Int64}
   use_dict_for_blocks::Bool
   fileindex_for_cellid_blocks::Dict{Int64,Int64} # [0] is index, [1] is blockcount
   cells_with_blocks::Dict{Int64,Int64}
   blocks_per_cell::Dict{Int64,Int64}
   blocks_per_cell_offsets::Dict{Int64,Int64}
   order_for_cellid_blocks::Dict{Int64,Int64}
   xcells::Int64
   ycells::Int64
   zcells::Int64
   xblock_size::Int64
   yblock_size::Int64
   zblock_size::Int64
   xmin::Float64
   ymin::Float64
   zmin::Float64
   xmax::Float64
   ymax::Float64
   zmax::Float64
   dx::Float64
   dy::Float64
   dz::Float64
   meshes::Dict{String, MeshInfo}
   populations::Vector{String}
end


function Base.show(io::IO, s::Meta)
   println(io, "filename = ", s.name)
end


"Return the xml footer of vlsv."
function read_xml_footer(fid)

   # First 8 bytes indicate big-endian or else
   endianness_offset = 8
   seek(fid, endianness_offset)
   # Obtain the offset of the XML file
   uint64_byte_amount = 8
   offset = read(fid, UInt64)
   seek(fid, offset)
   xmldata = read(fid, String)
   xmldoc  = parse_string(xmldata)
   footer = root(xmldoc)

end


"Return size and type information for the object."
function read_prep(fid, footer, name, tag, attr)

   arraysize = 0
   datasize = 0
   datatype = ""
   vectorsize = 0
   variable_offset = 0

   isFound = false
   
   for varinfo in footer[tag]
      if attribute(varinfo, attr) == name
         #has_attribute(e, name)
         arraysize = parse(Int, attribute(varinfo, "arraysize"))
         datasize = parse(Int, attribute(varinfo, "datasize"))
         datatype = attribute(varinfo, "datatype")
         vectorsize = parse(Int, attribute(varinfo, "vectorsize"))
         variable_offset = parse(Int, content(varinfo))
         isFound = true
         break
      end
   end

   if !isFound @error "unknown variable $(name)!" end

   seek(fid, variable_offset)

   if datatype == "float" && datasize == 4
      T = Float32
   elseif datatype == "float" && datasize == 8
      T = Float64
   elseif datatype == "int" && datasize == 4
      T = Int32
   elseif datatype == "int" && datasize == 8
      T = Int64
   elseif datatype == "uint" && datasize == 4
      T = UInt32
   elseif datatype == "uint" && datasize == 8
      T = UInt64
   end

   return T::DataType, variable_offset, arraysize, datasize, vectorsize
end


"""
Return the cells ID.
"""
read_fileindex_for_cellid(fid, footer) = 
   read_vector(fid, footer, "CellID", "VARIABLE")
   

"Return vector data from vlsv file."
function read_vector(fid, footer, name, tag)

   T, _, arraysize, _, _ = read_prep(fid, footer, name, tag, "name")

   w = Vector{T}(undef, arraysize)
   read!(fid, w)

   return w
end



function read_general(fid, footer; name="", tag="", mesh="", operator="pass", cellids=[])

end


function read_meta(filename::AbstractString; verbose=false)

   fid = open(filename, "r")

   use_dict_for_blocks = false
   fileindex_for_cellid_blocks = Dict()

   # Per population
   cells_with_blocks = Dict()
   blocks_per_cell = Dict()
   blocks_per_cell_offsets = Dict()
   order_for_cellid_blocks = Dict()

   footer = read_xml_footer(fid)

   meshName = "SpatialGrid"

   # [0] is index, [1] is blockcount
   fileindex_for_cellid = read_fileindex_for_cellid(fid, footer)

   cellIndex = sortperm(fileindex_for_cellid)

   bbox = read_mesh(fid, footer, meshName, "MESH_BBOX") 

   nodeCoordsX = read_mesh(fid, footer, meshName, "MESH_NODE_CRDS_X")
   nodeCoordsY = read_mesh(fid, footer, meshName, "MESH_NODE_CRDS_Y")
   nodeCoordsZ = read_mesh(fid, footer, meshName, "MESH_NODE_CRDS_Z")
  
   xcells, ycells, zcells = bbox[1:3]
   xblock_size, yblock_size, zblock_size = bbox[4:6]
   xmin, ymin, zmin = nodeCoordsX[1], nodeCoordsY[1], nodeCoordsZ[1]
   xmax, ymax, zmax = nodeCoordsX[end], nodeCoordsY[end], nodeCoordsZ[end]

   dx = (xmax - xmin) / xcells
   dy = (ymax - ymin) / ycells
   dz = (zmax - zmin) / zcells

   meshes = Dict{String,MeshInfo}()

   # Find all populations by the BLOCKIDS tag
   populations = String[]

   for varinfo in footer["BLOCKIDS"]

      if has_attribute(varinfo, "name")
         # New style vlsv file with bounding box
         popname = attribute(varinfo, "name")

         bbox = read_mesh(fid, footer, popname, "MESH_BBOX")

         nodeCoordsX = read_mesh(fid, footer, popname, "MESH_NODE_CRDS_X")   
         nodeCoordsY = read_mesh(fid, footer, popname, "MESH_NODE_CRDS_Y")   
         nodeCoordsZ = read_mesh(fid, footer, popname, "MESH_NODE_CRDS_Z")   
         vxblocks, vyblocks, vzblocks = bbox[1:3]
         vxblock_size, vyblock_size, vzblock_size = bbox[4:6]
         vxmin = nodeCoordsX[1]
         vymin = nodeCoordsY[1]
         vzmin = nodeCoordsZ[1]
         vxmax = nodeCoordsX[end]
         vymax = nodeCoordsY[end]
         vzmax = nodeCoordsZ[end]
         dvx = (vxmax - vxmin) / vxblocks / vxblock_size
         dvy = (vymax - vymin) / vyblocks / vyblock_size
         dvz = (vzmax - vzmin) / vzblocks / vzblock_size
      else
         popname = "avgs"

         if "vxblocks_ini" in attribute.(footer["PARAMETER"],"name") 
            # Old vlsv files where the mesh is defined with parameters
            vxblocks = read_parameter(fid, footer, "vxblocks_ini")
            vyblocks = read_parameter(fid, footer, "vyblocks_ini")
            vzblocks = read_parameter(fid, footer, "vzblocks_ini")
            vxblock_size = 4
            vyblock_size = 4
            vzblock_size = 4
            vxmin = read_parameter(fid, footer, "vxmin")
            vymin = read_parameter(fid, footer, "vymin")
            vzmin = read_parameter(fid, footer, "vzmin")
            vxmax = read_parameter(fid, footer, "vxmax")
            vymax = read_parameter(fid, footer, "vymax")
            vzmax = read_parameter(fid, footer, "vzmax")
            dvx = (vxmax - vxmin) / vxblocks / vxblock_size
            dvy = (vymax - vymin) / vyblocks / vyblock_size
            dvz = (vzmax - vzmin) / vzblocks / vzblock_size
         else
            # No velocity space info, e.g., file not written by Vlasiator 
            vxblocks, vyblocks, vzblocks = 0, 0, 0
            vxblock_size, vyblock_size, vzblock_size = 4, 4, 4
            vxmin, vymin, vzmin = 0.0, 0.0, 0.0
            vxmax, vymax, vzmax = 0.0, 0.0, 0.0
            dvx, dvy, dvz = 1.0, 1.0, 1.0
         end
      end

      # Update list of active populations
      if popname ∉ populations 
         push!(populations, popname)
      end

      # Create a new MeshInfo object for this population
      popMesh = MeshInfo(vxblocks, vyblocks, vzblocks, 
         vxblock_size, vyblock_size, vzblock_size,
         vxmin, vymin, vzmin, vxmax, vymax, vzmax,
         dvx, dvy, dvz)

      meshes[popname] = popMesh

      if verbose
         @info "Found population " * popname
      end
   end

   #close(fid) # Is it safe not to close it?

   meta = Meta(filename, fid, footer, fileindex_for_cellid, cellIndex,
      use_dict_for_blocks, 
      fileindex_for_cellid_blocks, cells_with_blocks, blocks_per_cell, 
      blocks_per_cell_offsets, order_for_cellid_blocks, 
      xcells, ycells, zcells, xblock_size, yblock_size, zblock_size,
      xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, meshes, 
      populations)

end


"""
   read_variable_info(footer, var)

Return a struct of VarInfo.           
"""
function read_variable_info(footer, var)

   unit = ""
   unitLaTeX = ""
   variableLaTeX = ""
   unitConversion = ""

   # Force lowercase
   var = lowercase(var)

   # Get population and variable names from data array name 
   if occursin("/", var)
      popname, varname = split(a, "/")
   else
      popname = "pop"
      varname = var
   end
      
#   # Check which set of datareducers to use
#   if varname[1:3] == "vg_"
#      reducer_reg = v5reducers
#      reducer_multipop = multipopv5reducers
#   else
#      reducer_reg = datareducers
#      reducer_multipop = multipopdatareducers
#   end

   if has_variable(footer, var)
      if varname[1:3] == "vg_" || varname[1:3] == "fg_"
         # For Vlasiator 5 vlsv files, metadata is included

         for varinfo in footer["VARIABLE"]
            if attribute(varinfo, "name") == var
               unit = attribute(varinfo, "unit")
               unitLaTeX = attribute(varinfo, "unitLaTeX")
               variableLaTeX = attribute(varinfo, "variableLaTeX")
               unitConversion = attribute(varinfo, "unitConversion") 
            end
         end

         # Correction for early version incorrect number density (extra backslash)
         if variableLaTeX[1:3] == r"$\n"
            variableLaTeX = r"$n"*variableLaTeX[4:end]
         end

      elseif var in keys(units_predefined)
         unit = units_predefined[var]
         variableLaTeX = latex_predefined[var]
         unitLaTeX = latexunits_predefined[var]
      end    
#   elseif var in reducer_reg
#      unit = reducer_reg[name].unit
#      varLaTeX = reducer_reg[name].varLaTeX
#      unitLaTeX = reducer_reg[name].unitLaTeX
#   elseif "pop/"*varname in reducer_multipop
#      poplatex = "i"
#      if popname in vlsvvariables.speciesdict
#         poplatex = vlsvvariables.speciesdict[popname]
#      end
#      unit = reducer_multipop["pop/"*varname].unit
#      varLaTeX = (reducer_multipop["pop/"*varname].varLaTeX).replace("REPLACEPOP",poplatex)
#      unitLaTeX = reducer_multipop["pop/"*varname].unitLaTeX
#   elseif varname in vlsvvariables.unitsdict
#      poplatex = "i"
#      if popname in vlsvvariables.speciesdict
#         poplatex = vlsvvariables.speciesdict[check_variablepopname]
#      end
#      unit = vlsvvariables.unitsdict[varname]
#      varLaTeX = vlsvvariables.latexdictmultipop[varname].replace("REPLACEPOP",poplatex)
#      unitLaTeX = vlsvvariables.latexunitsdict[varname]
   else
      #unit = ""
      #varLaTeX = r""*name.replace("_","\_")
      #unitLaTeX = ""
   end

#   if startswith(name, "fg_")
#       data = read_fsgrid_variable(name=name, operator=operator)
#   else
#       data = read_variable(name=name, operator=operator, cellids=cellids)
#   end

#   if operator != "pass"
#      if operator == "magnitude"
#         varLaTeX = r"$|$"*varLaTeX*r"$|$"
#      else
#         varLaTeX = varLaTeX*r"${_{"*operator*r"}}$"
#      end
#
#   end

   return VarInfo(unit, unitLaTeX, variableLaTeX, unitConversion)
end


function read_mesh(fid, footer, typeMesh, varMesh)

   T, variable_offset, arraysize, datasize, vectorsize = 
      read_prep(fid, footer, typeMesh, varMesh, "mesh")

   w = Vector{T}(undef, arraysize)
   read!(fid, w)

   return w
end

"""
   read_variable(meta, var)

Return variable var from the vlsv file.
"""
read_variable(meta, var) = 
   read_vector(meta.fid, meta.footer, var, "VARIABLE")


has_variable(footer, var) = has_name(footer, "VARIABLE", var)


function read_variable_select(meta, var, cellIDs=UInt[])

   if isempty(cellIDs)
      w = read_variable(meta, var)
      return [w]
   end

   T, variable_offset, arraysize, datasize, vectorsize = 
      read_prep(meta.fid, meta.footer, var, "VARIABLE", "name")

   rOffsets = [meta.cellIndex[i]*datasize*vectorsize for i in cellIDs]

   v = fill(T[], length(rOffsets))

   for (i, r) in enumerate(rOffsets)
      loc = variable_offset + r
      seek(meta.fid, loc)
   
      w = Vector{T}(undef, vectorsize)
      read!(meta.fid, w)
      v[i] = w
   end

   return v
end


function read_parameter(fid, footer, param)
   
   T, _, _, _, _ = read_prep(fid, footer, param, "PARAMETER", "name")

   p = read(fid, T)
end


read_parameter(meta, param) = read_parameter(meta.fid, meta.footer, param)


has_parameter(meta, param) = has_name(meta.footer, "PARAMETER", param)


function has_name(footer, tag, name)
   isFound = false
   
   for varinfo in footer[tag]
      attribute(varinfo, "name") == name && (isFound = true)
      isFound && break
   end
   
   return isFound
end