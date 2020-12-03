# Vlasiator plotting in Julia.
#
# Hongyang Zhou, hyzhou@umich.edu 12/03/2020

using PyPlot

function plot_contour(meta, var)

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
   elseif ysize == 1 && zsize ==1

   end

   datainfo = read_var_meta(meta.footer, var)

   data = read_variable(meta, var)

   # fsgrid array
   if startswith(var, "fg_")
      #data = swapaxes(data, 0,1)
   else # vlasov grid
      if ndims(data) == 1          
          data = reshape(data[meta.cellIndex], sizes[1], sizes[2])
      elseif ndims(data) == 2
          #data = data[cellids.argsort()].reshape([sizes[2],sizes[1],data.shape[1]])
      elseif ndims(data) == 3
          #data = data[cellids.argsort()].reshape([sizes[2],sizes[1],data.shape[1],data.shape[2]])
      else
          @error "Error in reshaping data $(var)!"
      end
   end

   x = range(plotrange[1], plotrange[2], length=sizes[1])
   y = range(plotrange[1], plotrange[2], length=sizes[2])

   # Logarithmic plot
   vmin = minimum(data[data .> 0.0])
   vmax = maximum(data)

   norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)
   ticks = matplotlib.ticker.LogLocator(base=10,subs=collect(0:9))
   cmap = matplotlib.cm.turbo

   fig, ax = subplots()

   c = ax.pcolormesh(x, y, data', norm=norm, cmap=cmap, shading="auto")

   fig.colorbar(c, ticks=ticks)

end
