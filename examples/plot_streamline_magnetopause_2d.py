"""
An example of using find_magnetopause_2d() and plotting the output over colormap.
"""

import analysator as pt

run = 'BFD'
file_name="/wrk-vakka/group/spacephysics/vlasiator/2D/"+run+"/bulk/bulk.0002000.vlsv"
outdir= ""

# magnetopause points
magnetopause = pt.calculations.find_magnetopause_sw_streamline_2d(
        file_name,
        streamline_seeds=None,
        seeds_n=200,
        seeds_x0=20*6371000,
        seeds_range=[-5*6371000, 5*6371000],
        streamline_length=45*6371000,
        end_x=-15*6371000,
        x_point_n=50)

magnetopause= magnetopause*(1/6371000) # in RE

def external_plot(ax,XmeshXY=None, YmeshXY=None, pass_maps=None):
        ax.plot(magnetopause[:,0], magnetopause[:,1], color='cyan', linewidth=1.0) 

pt.plot.plot_colormap(
    filename=file_name,
    var="rho",
    boxre = [-21,21,-21,21],
    Earth=True,
    external = external_plot,
    run=run,
    colormap='inferno',
    )
