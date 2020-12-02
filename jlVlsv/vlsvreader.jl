# VLSV reader in Julia
#
# Hongyang Zhou, hyzhou@umich.edu 12/01/2020

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


"Meta data declaration."
struct Meta
   name::AbstractString
   fid::IOStream
   footer::XMLElement
   fileindex_for_cellid::Dict{UInt64,Int64}
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
   
   for varinfo in footer[tag]
      if attribute(varinfo, attr) == name
         arraysize = parse(Int, attribute(varinfo, "arraysize"))
         datasize = parse(Int, attribute(varinfo, "datasize"))
         datatype = attribute(varinfo, "datatype")
         vectorsize = parse(Int, attribute(varinfo, "vectorsize"))
         variable_offset = parse(Int, content(varinfo))
         break
      end
   end

   if variable_offset == 0
      @error "unknown variable $(name)!"
   end

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
Read in the cell ids and create an internal dictionary to give the index of an
arbitrary cellID.
"""
function read_fileindex_for_cellid(fid, footer)

   w = read_vector(fid, footer, "CellID", "VARIABLE")

   fileindex_for_cellid = 
      Dict(cellid => index for (index,cellid) in enumerate(w))
   
end

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
      popname = attribute(varinfo, "name"; required=false)
      if isnothing(popname) 
         popname = "avgs"
      end

      # Update list of active populations
      if popname ∉ populations 
         push!(populations, popname)
      end

      bbox = read_mesh(fid, footer, popname, "MESH_BBOX")

      if !isempty(bbox) # new style vlsv file with bounding box
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
         # Velocity cell lengths
         dvx = (vxmax - vxmin) / vxblocks / vxblock_size
         dvy = (vymax - vymin) / vyblocks / vyblock_size
         dvz = (vzmax - vzmin) / vzblocks / vzblock_size
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

   meta = Meta(filename, fid, footer, fileindex_for_cellid, use_dict_for_blocks, 
      fileindex_for_cellid_blocks, cells_with_blocks, blocks_per_cell, 
      blocks_per_cell_offsets, order_for_cellid_blocks, 
      xcells, ycells, zcells, xblock_size, yblock_size, zblock_size,
      xmin, ymin, zmin, xmax, ymax, zmax, dx, dy, dz, meshes, 
      populations)

end


"""
   read_var_meta(footer, var)

Return variable metadata from the vlsv footer, including the unit of the 
variable as a regular string, the unit of the variable as a LaTeX-formatted 
string, the description of the variable as a LaTeX-formatted string, and the 
conversion factor to SI units as a string.             
"""
function read_var_meta(footer, var)

   unit = ""
   unitLaTeX = ""
   variableLaTeX = ""
   unitConversion = ""

   for varinfo in footer["VARIABLE"]
      if attribute(varinfo, "name") == var
         unit = attribute(varinfo, "unit")
         unitLaTeX = attribute(varinfo, "unitLaTeX")
         variableLaTeX = attribute(varinfo, "variableLaTeX")
         unitConversion = attribute(varinfo, "unitConversion") 
      end
   end

   return unit, unitLaTeX, variableLaTeX, unitConversion
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


function read_variable_select(meta, var, cellIDs=UInt[])

   if isempty(cellIDs)
      w = read_variable(meta, var)
      return [w]
   end

   T, variable_offset, arraysize, datasize, vectorsize = 
      read_prep(meta.fid, meta.footer, var, "VARIABLE", "name")

   rOffsets = [meta.fileindex_for_cellid[i]*datasize*vectorsize for i in cellIDs]

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


function read_parameter(meta, param)

   fid, footer = meta.fid, meta.footer

   T, _, _, _, _ = read_prep(fid, footer, param, "PARAMETER", "name")

   p = read(fid, T)

end
