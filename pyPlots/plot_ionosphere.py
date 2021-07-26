import matplotlib
import pytools as pt
import numpy as np
import matplotlib.pyplot as plt
import glob
import os, sys
import time
import re
import matplotlib.ticker as mtick
from matplotlib.ticker import MaxNLocator, MultipleLocator
from matplotlib.ticker import LogLocator
from matplotlib.colors import BoundaryNorm,LogNorm,SymLogNorm
from distutils.version import LooseVersion, StrictVersion

def plot_ionosphere(filename=None,
                  vlsvobj=None,
                  filedir=None, step=None, run=None,
                  outputdir=None, outputfile=None,
                  nooverwrite=False,
                  var=None, op=None, operator=None,
                  colormap=None, vmin=None, vmax=None,
                  symmetric=False, absolute=None,
                  usesci=True,
                  lin=True, symlog=None, nocb=False,
                  minlatitude=50,
                  cbtitle=None, title=None,
                  thick=1.0,scale=1.0,vscale=1.0,
                  viewdir=1.0,draw=None
      ):

    ''' Plots a 2d projection of an ionosphere variable

    :kword filename:    path to .vlsv file to use for input. Assumes a bulk file.
    :kword vlsvobj:     Optionally provide a python vlsvfile object instead
    :kword filedir:     Optionally provide directory where files are located and use step for bulk file name
    :kword step:        output step index, used for constructing output (and possibly input) filename
    :kword run:         run identifier, used for constructing output filename
    :kword outputdir:   path to directory where output files are created (default: $HOME/Plots/ or override with PTOUTPUTDIR)
                        If directory does not exist, it will be created. If the string does not end in a
                        forward slash, the final part will be used as a prefix for the files.
    :kword outputfile:  Singular output file name
    :kword nooverwrite: Set to only perform actions if the target output file does not yet exist                    
    :kword var:         variable to plot, e.g. rho, RhoBackstream, beta, Temperature, MA, Mms, va, vms,
                        E, B, v, V or others. Accepts any variable known by analysator/pytools.
                        Per-population variables are simply given as "proton/rho" etc
    :kword operator:    Operator to apply to variable: None, x, y, or z. Vector variables return either
                        the queried component, or otherwise the magnitude. 
    :kword op:          duplicate of operator
    :kword colormap:    colour scale for plot, use e.g. hot_desaturated, jet, viridis, plasma, inferno,
                        magma, parula, nipy_spectral, RdBu, bwr
    :kword vmin,vmax:   min and max values for colour scale and colour bar. If no values are given,
                        min and max values for whole plot (non-zero rho regions only) are used.
    :kword symmetric:   Set the absolute value of vmin and vmax to the greater of the two
    :kword absolute:    Plot the absolute of the evaluated variable
    :kwird usesci:      Use scientific notation for colorbar ticks? (default: True)
    :kword lin:         Flag for using linear colour scaling instead of log
    :kword symlog:      Use logarithmic scaling, but linear when abs(value) is below the value given to symlog.
                        Allows symmetric quasi-logarithmic plots of e.g. transverse field components.
                        A given of 0 translates to a threshold of max(abs(vmin),abs(vmax)) * 1.e-2, but this can
                        result in the innermost tick marks overlapping. In this case, using a larger value for 
                        symlog is suggested.
    :kword nocb:        Set to suppress drawing of colourbar
    :kword minlatitude: Minimum plot latitude (default=50 degrees)
    :kword title:       string to use as plot title instead of time.
                        Special case: Set to "msec" to plot time with millisecond accuracy or "musec"
                        for microsecond accuracy. "sec" is integer second accuracy.
    :kword cbtitle:     string to use as colorbar title instead of map name
    :kword viewdir:     view direction onto the sphere. Positive value: North pole. Negative values: South pole.
    :kword scale:       Scale text size (default=1.0)
    :kword vscale:      Scale all values with this before plotting. Useful for going from e.g. m^-3 to cm^-3
                        or from tesla to nanotesla. Guesses correct units for colourbar for some known
                        variables.
    :kword thick:       line and axis thickness, default=1.0
    :kword draw:        Set to anything but None or False in order to draw image on-screen instead of saving to file (requires x-windowing)

 
    .. code-block:: python

    '''

    IONOSPHERE_RADIUS=6471e3 # R_E + 100 km
    # TODO: Implement vlasiator logo watermarking

    # Change certain falsy values:
    if not lin and lin != 0:
        lin = None
    if not symlog and symlog != 0:
        symlog = None
    if symlog is True:
        symlog = 0
    if (filedir == ''):
        filedir = './'
    if (outputdir == ''):
        outputdir = './'

    # Input file or object
    if filename!=None:
        f=pt.vlsvfile.VlsvReader(filename)
    elif ((filedir!=None) and (step!=None)):
        filename = glob.glob(filedir+'bulk*'+str(step).rjust(7,'0')+'.vlsv')[0]
        f=pt.vlsvfile.VlsvReader(filename)
    elif vlsvobj!=None:
        f=vlsvobj
    else:
        print("Error, needs a .vlsv file name, python object, or directory and step")
        return

    if operator is None:
        if op is not None:
            operator=op

    if not colormap:
        # Default values
        colormap="bwr"
    cmapuse=matplotlib.cm.get_cmap(name=colormap)

    fontsize=8*scale # Most text
    fontsize2=10*scale # Time title
    fontsize3=8*scale # Colour bar ticks and title

    # Plot title with time
    timeval=f.read_parameter("time")

    # Plot title with time
    if title is None or title=="msec" or title=="musec":        
        if timeval is None:    
            plot_title = ''
        else:
            timeformat='{:4.1f}'
            if title=="sec": timeformat='{:4.0f}'
            if title=="msec": timeformat='{:4.3f}'
            if title=="musec": timeformat='{:4.6f}'
            plot_title = "t="+timeformat.format(timeval)+' s '
    else:
        plot_title = ''+title

    # step, used for file name
    if step is not None:
        stepstr = '_'+str(step).rjust(7,'0')
    else:
        if filename:
            stepstr = '_'+filename[-12:-5]
        else:
            stepstr = ''

    # If run name isn't given, just put "plot" in the output file name
    if run is None:
        run='plot'
        if filename:
            # If working within CSC filesystem, make a guess:
            if filename[0:16]=="/proj/vlasov/3D/":
                run = filename[16:19]

    # Verify validity of operator
    operatorstr=''
    operatorfilestr=''
    if operator is not None:
        # .isdigit checks if the operator is an integer (for taking an element from a vector)
        if type(operator) is int:
            operator = str(operator)
        if not operator in 'xyz' and operator != 'magnitude' and not operator.isdigit():
            print("Unknown operator "+str(operator))
            operator=None
            operatorstr=''
        if operator in 'xyz':
            # For components, always use linear scale, unless symlog is set
            operatorstr='_'+operator
            operatorfilestr='_'+operator
            if symlog is None:
                lin=True
        # index a vector
        if operator.isdigit():
            operator = str(operator)
            operatorstr='_{'+operator+'}'
            operatorfilestr='_'+operator
        # Note: operator magnitude gets operatorstr=''

    # Output file name
    if not var:
        # If no expression or variable given, defaults to FACs
        var='ig_fac'
    varstr=var.replace("/","_")

    if not outputfile: # Generate filename
        if not outputdir: # default initial path
            outputdir=pt.plot.defaultoutputdir
        # Sub-directories can still be defined in the "run" variable
        outputfile = outputdir+run+"_ionosphere_"+varstr+operatorfilestr+stepstr+".png"
    else: 
        if outputdir:
            outputfile = outputdir+outputfile

    # Re-check to find actual target sub-directory
    outputprefixind = outputfile.rfind('/')
    if outputprefixind >= 0:            
        outputdir = outputfile[:outputprefixind+1]

    # Ensure output directory exists
    if not os.path.exists(outputdir):
        try:
            os.makedirs(outputdir)
        except:
            pass

    if not os.access(outputdir, os.W_OK):
        print(("No write access for directory "+outputdir+"! Exiting."))
        return

    # Check if target file already exists and overwriting is disabled
    if (nooverwrite and os.path.exists(outputfile)):            
        if os.stat(outputfile).st_size > 0: # Also check that file is not empty
            print(("Found existing file "+outputfile+". Skipping."))
            return
        else:
            print(("Found existing file "+outputfile+" of size zero. Re-rendering."))

    # The plot will be saved in a new figure ('draw' and 'axes' keywords not implemented yet)
    if str(matplotlib.get_backend()) != pt.backend_noninteractive: #'Agg':
        plt.switch_backend(pt.backend_noninteractive)

    # Select plotting back-end based on on-screen plotting or direct to file without requiring x-windowing
    if draw:
        if str(matplotlib.get_backend()) != pt.backend_interactive: #'TkAgg': 
            plt.switch_backend(pt.backend_interactive)
    else:
        if str(matplotlib.get_backend()) != pt.backend_noninteractive: #'Agg':
            plt.switch_backend(pt.backend_noninteractive)  

    # Read ionosphere mesh node coordinates
    coords = f.get_ionosphere_node_coords()
    # Read ionosphere mesh connectivity
    elements = f.get_ionosphere_element_corners()
    # Read selected variable valus
    if operator is None:
       operator="pass"
    datamap_info= f.read_variable_info(var, operator=operator)
    # Correction for broken ionosphere units
    datamap_info.latex = re.sub("\\\\text","\\\\mathrm", datamap_info.latex)
    datamap_info.latexunits = re.sub("\\\\mho","\\\\Omega^{-1}", datamap_info.latexunits)
    cb_title_use = datamap_info.latex
    datamap_unit = datamap_info.latexunits
    # Check if vscale results in standard unit
    datamap_unit = pt.plot.scaleunits(datamap_info, vscale)
    values = datamap_info.data

    # Add unit to colorbar title
    if datamap_unit:
       cb_title_use = cb_title_use + "\,["+datamap_unit+"]"

    # Project nodes and elements into view plane
    mask = []
    if viewdir > 0:
      for e in elements:
         if coords[e[0],2] > 0 and coords[e[1],2] > 0 and coords[e[2],2] > 0:
            mask+=[False]
         else:
            mask+=[True]
    else:
      for e in elements:
         if coords[e[0],2] < 0 and coords[e[1],2] < 0 and coords[e[2],2] < 0:
            mask+=[False]
         else:
            mask+=[True]

    # Build mesh triangulation
    r=np.degrees(np.arccos(coords[:,2]/IONOSPHERE_RADIUS))
    theta=np.arctan2(coords[:,1],coords[:,0])
    tri = matplotlib.tri.Triangulation(r*np.cos(theta), r*np.sin(theta), elements, mask)

    # Allow title override
    if cbtitle is not None:
        # Here allow underscores for manual math mode
        cb_title_use = cbtitle

    # If automatic range finding is required, find min and max of array
    if vmin is not None:
        vminuse=vmin
    else: 
        vminuse = np.ma.amin(values)
    if vmax is not None:
        vmaxuse=vmax
    else:
        vmaxuse = np.ma.amax(values)

    # If both values are zero, we have an empty array
    if vmaxuse==vminuse==0:
        print("Error, requested array is zero everywhere. Exiting.")
        return 0

    # If vminuse and vmaxuse are extracted from data, different signs, and close to each other, adjust to be symmetric
    # e.g. to plot transverse field components. Always done for symlog.
    if vmin is None and vmax is None:
        if symmetric or np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2) or symlog is not None:
            absval = max(abs(vminuse),abs(vmaxuse))
            vminuse = -absval
            vmaxuse = absval

    # Ensure that lower bound is valid for logarithmic plots
    if (vminuse <= 0) and (lin is None) and (symlog is None):
        # Drop negative and zero values
        vminuse = min(np.ma.amin(datamap), 0)

    # Special case of very small vminuse values
    if ((vmin is None) or (vmax is None)) and (vminuse > 0) and (vminuse < vmaxuse*1.e-5):
        vminuse = vmaxuse*1e-5
        if lin is not None:
            vminuse = 0

    # If symlog scaling is set:
    if symlog is not None:
        if symlog > 0:
            linthresh = symlog 
        else:
            linthresh = max(abs(vminuse),abs(vmaxuse))*1.e-2

    # Lin or log colour scaling, defaults to log
    if lin is None:
        # Special SymLogNorm case
        if symlog is not None:
            if LooseVersion(matplotlib.__version__) < LooseVersion("3.2.0"):
                norm = SymLogNorm(linthresh=linthresh, linscale = 1.0, vmin=vminuse, vmax=vmaxuse, clip=True)
                print("WARNING: colormap SymLogNorm uses base-e but ticks are calculated with base-10.")
                #TODO: copy over matplotlib 3.3.0 implementation of SymLogNorm into pytools/analysator
            else:
                norm = SymLogNorm(base=10, linthresh=linthresh, linscale = 1.0, vmin=vminuse, vmax=vmaxuse, clip=True)
            maxlog=int(np.ceil(np.log10(vmaxuse)))
            minlog=int(np.ceil(np.log10(-vminuse)))
            logthresh=int(np.floor(np.log10(linthresh)))
            logstep=1
            ticks=([-(10**x) for x in range(logthresh, minlog+1, logstep)][::-1]
                    +[0.0]
                    +[(10**x) for x in range(logthresh, maxlog+1, logstep)] )
        else:
            # Logarithmic plot
            norm = LogNorm(vmin=vminuse,vmax=vmaxuse)
            ticks = LogLocator(base=10,subs=range(10)) # where to show labels
    else:
        # Linear
        levels = MaxNLocator(nbins=255).tick_values(vminuse,vmaxuse)
        norm = BoundaryNorm(levels, ncolors=cmapuse.N, clip=True)
        ticks = np.linspace(vminuse,vmaxuse,num=11)

    # Creating a new figure and axes 
    figsize = (6,5)
    if nocb:
        figsize = (5,5)
    fig = plt.figure(figsize=figsize,dpi=150)
    ax_cartesian = fig.add_axes([0.1,0.1,0.9,0.9], xlim=(-(90-minlatitude),(90-minlatitude)), ylim=(-(90-minlatitude),(90-minlatitude)), aspect='equal')
    ax_cartesian.set_xticklabels([])
    ax_cartesian.set_yticklabels([])
    ax_cartesian.axis('off')

    ax_polar = fig.add_axes([0.1,0.1,0.9,0.9], polar=True, frameon=False, ylim=(0, minlatitude))


    ### THE ACTUAL PLOT HAPPENS HERE ###
    contours = ax_cartesian.tricontourf(tri, values, cmap=cmapuse, vmin=vmin, vmax=vmax, levels=64)
    ax_polar.grid(True)
    ax_polar.set_rgrids(range(0,minlatitude,10), map(lambda x: str(90-x)+"°", range(0,minlatitude,10)),angle=310)


    # Colourbar title
    if len(cb_title_use)!=0:
        cb_title_use = pt.plot.mathmode(pt.plot.bfstring(cb_title_use))

    # Set flag which affects colorbar decimal precision
    if lin is None:
        pt.plot.cb_linear = False
    else:
        pt.plot.cb_linear = True

    # Creating colorbar axes
    if not nocb:
        cax = fig.add_axes([0.76,0.2,0.03,0.6])
        cbdir="right"; horalign="left"

        # First draw colorbar
        if usesci:
            cb = plt.colorbar(contours, ticks=ticks, format=mtick.FuncFormatter(pt.plot.cbfmtsci), cax=cax, drawedges=False)
        else:
            cb = plt.colorbar(contours, ticks=ticks, format=mtick.FuncFormatter(pt.plot.cbfmt), cax=cax, drawedges=False)
        cb.outline.set_linewidth(thick)
        cb.ax.yaxis.set_ticks_position(cbdir)
 
        cbticks = cb.get_ticks()
        cb.set_ticks(cbticks[(cbticks>=vminuse)*(cbticks<=vmaxuse)])

        cb.ax.tick_params(labelsize=fontsize3, width=thick)#,width=1.5,length=3)
        cb_title = cax.set_title(cb_title_use,fontsize=fontsize3,fontweight='bold', horizontalalignment=horalign)
        cb_title.set_position((0.,1.+0.025*scale)) # avoids having colourbar title too low when fontsize is increased

        # Perform intermediate draw if necessary to gain access to ticks
        if (symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2)) or (not lin and symlog is None):
            fig.canvas.draw() # draw to get tick positions

        # Adjust placement of innermost ticks for symlog if it indeed is (quasi)symmetric
        if symlog is not None and np.isclose(vminuse/vmaxuse, -1.0, rtol=0.2):
            cbt=cb.ax.yaxis.get_ticklabels()
            (cbtx,cbty) = cbt[len(cbt)//2-1].get_position() # just below zero
            if abs(0.5-cbty)/scale < 0.1:
                cbt[len(cbt)//2-1].set_va("top")
            (cbtx,cbty) = cbt[len(cbt)//2+1].get_position() # just above zero
            if abs(0.5-cbty)/scale < 0.1:
                cbt[len(cbt)//2+1].set_va("bottom")
            if len(cbt)>=7: # If we have at least seven ticks, may want to adjust next ones as well
                (cbtx,cbty) = cbt[len(cbt)//2-2].get_position() # second below zero
                if abs(0.5-cbty)/scale < 0.15:
                    cbt[len(cbt)//2-2].set_va("top")
                (cbtx,cbty) = cbt[len(cbt)//2+2].get_position() # second above zero
                if abs(0.5-cbty)/scale < 0.15:
                    cbt[len(cbt)//2+2].set_va("bottom")

        # Adjust precision for colorbar ticks
        thesetickvalues = cb.locator()
        if len(thesetickvalues)<2:
            precision_b=1
        else:
            mintickinterval = abs(thesetickvalues[-1]-thesetickvalues[0])
            # find smallest interval
            for ticki in range(len(thesetickvalues)-1):
                mintickinterval = min(mintickinterval,abs(thesetickvalues[ticki+1]-thesetickvalues[ticki]))
            precision_a, precision_b = '{:.1e}'.format(mintickinterval).split('e')
            # e.g. 9.0e-1 means we need precision 1
            # e.g. 1.33e-1 means we need precision 3?
        pt.plot.decimalprecision_cblin = 1
        if int(precision_b)<1: pt.plot.decimalprecision_cblin = str(1+abs(-int(precision_b)))
        cb.update_ticks()

        # if too many subticks in logarithmic colorbar:
        if not lin and symlog is None:
            nlabels = len(cb.ax.yaxis.get_ticklabels()) # TODO implement some kind of ratio like in other scripts, if needed?
            valids = ['1','2','3','4','5','6','7','8','9']
            if nlabels > 10:
                valids = ['1','2','3','4','5','6','8']
            if nlabels > 19:
                valids = ['1','2','5']
            if nlabels > 28:
                valids = ['1']
            # for label in cb.ax.yaxis.get_ticklabels()[::labelincrement]:
            for labi,label in enumerate(cb.ax.yaxis.get_ticklabels()):
                labeltext = label.get_text().replace('$','').replace('{','').replace('}','').replace('\mbox{\textbf{--}}','').replace('-','').replace('.','').lstrip('0')
                if not labeltext:
                    continue
                firstdigit = labeltext[0]
                if not firstdigit in valids: 
                    label.set_visible(False)

    plt.tight_layout(pad=0.01) 
    savefig_pad=0.01
    bbox_inches='tight'

    # Save output or draw on-screen
    if not draw:
        try:
            plt.savefig(outputfile,dpi=300, bbox_inches=bbox_inches, pad_inches=savefig_pad)
        except:
            print("Error attempting to save figure: ", sys.exc_info())
        print('...Done!') 
    else:
        # Draw on-screen
        plt.draw()
        plt.show()
        print('Draw complete!')
