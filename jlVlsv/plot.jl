# Vlasiator plotting in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu 12/03/2020

using PyPlot, Printf

const Re = 6.371e6 # Earth radius, [m]


"""
   plot_pcolormesh(meta, var (...))

Plot a 2D pseudocolor var from vlsv.

`plot_pcolormesh(meta, var, )`

`plot_pcolormesh(meta, var, axisunit="Re)`

`plot_pcolormesh(data, func, isLinear=false)`
"""
function plot_pcolormesh(meta, var; axisunit="Re", isLinear=false)

   xsize, ysize, zsize = meta.xcells, meta.ycells, meta.zcells
   xmin, ymin, zmin = meta.xmin, meta.ymin, meta.zmin 
   xmax, ymax, zmax = meta.xmax, meta.ymax, meta.zmax
   dx = meta.dx # cell size is equal in x,y,z for now

   # Check if ecliptic or polar run
   if ysize == 1 && zsize != 1
      plotrange = [xmin,xmax,zmin,zmax]
      sizes = [xsize,zsize]
      PLANE = "XZ"
   elseif zsize == 1 && ysize != 1
      plotrange = [xmin,xmax,ymin,ymax]
      sizes = [xsize,ysize]
      PLANE = "XY"
   elseif ysize == 1 && zsize == 1

   end

   datainfo = read_variable_info(meta.footer, var)

   if var in keys(variables_predefined)
      data = variables_predefined[var](meta)
   else
      data = read_variable(meta, var)
   end

   # fsgrid array
   if startswith(var, "fg_")
      #data = swapaxes(data, 0, 1)
   else # vlasov grid
      if ndims(data) == 1 || (ndims(data) == 2 && size(data)[1] == 1)       
         data = reshape(data[meta.cellIndex], sizes[1], sizes[2])
      elseif ndims(data) == 2
         #data = data[cellids.argsort()].reshape([sizes[2],sizes[1],data.shape[1]])
      elseif ndims(data) == 3
         #data = data[cellids.argsort()].reshape([sizes[2],sizes[1],data.shape[1],data.shape[2]])
      else
         @error "Error in reshaping data $(var)!"
      end
   end

   if axisunit == "Re"
      x = range(plotrange[1], plotrange[2], length=sizes[1]) ./ Re
      y = range(plotrange[1], plotrange[2], length=sizes[2]) ./ Re      
   else
      x = range(plotrange[1], plotrange[2], length=sizes[1])
      y = range(plotrange[1], plotrange[2], length=sizes[2])
   end

   if any([has_parameter(meta, p) for p in ("t", "time")])
      timesim = read_parameter(meta, "t")   
      plot_title = @sprintf "t= %4.1fs" timesim
   else
      plot_title = ""
   end

   cmap = matplotlib.cm.turbo

   if !isLinear
      # Logarithmic plot
      vmin = minimum(data[data .> 0.0])
      vmax = maximum(data)
      
      norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
      ticks = matplotlib.ticker.LogLocator(base=10,subs=collect(0:9))
   else
      vmin = minimum(data)
      vmax = maximum(data)
      nticks = 7
      levels = matplotlib.ticker.MaxNLocator(nbins=255).tick_values(vmin, vmax)
      norm = matplotlib.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=true)
      ticks = range(vmin, vmax, length=nticks)
   end

   cb_title_use = datainfo.variableLaTeX
   data_unit = datainfo.unitLaTeX
   cb_title_use *= ",["*data_unit*"]"

   fig, ax = subplots()

   c = ax.pcolormesh(x, y, data', norm=norm, cmap=cmap, shading="auto")

   cb = fig.colorbar(c, ticks=ticks)

   ax.set_title(plot_title, fontsize=14, fontweight="bold")
   
   ax.set_aspect("equal")

   for axis in ["top","bottom","left","right"]
       ax.spines[axis].set_linewidth(2.0)
   end
   ax.xaxis.set_tick_params(width=2.0,length=3)
   ax.yaxis.set_tick_params(width=2.0,length=3)
   
   cb_title = cb.ax.set_title(cb_title_use, fontsize=14, fontweight="bold")
   cb.outline.set_linewidth(2.0)

   return c
end