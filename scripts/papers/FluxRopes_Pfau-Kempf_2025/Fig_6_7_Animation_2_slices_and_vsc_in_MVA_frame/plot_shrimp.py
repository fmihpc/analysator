#!/usr/bin/env python
import matplotlib
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import pylab as pl
import numpy as np
import math
import sys
import pytools as pt
from variable import get_data, get_name, get_units


#plt.xkcd(scale=1.5, length=100, randomness=3)

plot_lmn = True
circleradius = 0 # 35.5
circlex = 18
circley = 0

#import plot_colormap3dslice_withcontours

contourvariable = "vg_fluxrope"
contouroperator = "pass"
contourvalues = 7
contourcolors = "xkcd:swamp"
contourlinewidths = 7

lw=4 # VSC plots
thick=3
panel_labels = ("(c)", "(d)", "(e)", "(f)", "(a)", "(b)")

Re = 6.371e+6 # Earth radius in m

numVSC = 4
colorseq = ("xkcd:gunmetal", "cyan", "brown", "xkcd:radioactive green")
bgrseq = ('blue', '#E69F00', 'red')
lmn = ('l', 'm', 'n')
linestyleseq = ('dashed', 'dashdot', 'dotted')
# VSC positions in metres
pos = []

pos.append([-6.4e+07, -17*Re, -3*Re])
pos.append([-6.4e+07, -17*Re, -4.5*Re])
pos.append([-6.4e+07, -17*Re, -6*Re])
pos.append([-6.4e+07, -17*Re, -7.5*Re])

pos = np.array(pos)

angle = np.arctan2(pos[:,2]/Re - circley, pos[:,1]/Re - circlex)

arrowwidth = 0.01
headwidth_default = 3
headlength_default = 5
headaxislength_default = 5
triad_zorder = 1e6

# initialize lists for gathering values
profiles_B = []
profiles_Blmnmag = []
profiles_Bx = []
profiles_By = []
profiles_Bz = []
profiles_Bl = []
profiles_Bm = []
profiles_Bn = []
profiles_V = []
profiles_Vx = []
profiles_Vy = []
profiles_Vz = []
profiles_rho = []
profiles_T = []
profiles_times = []
vsc_eigenvalues = []
vsc_eigenvectors = []

for i in range(numVSC):
    profiles_B.append([])
    profiles_Blmnmag.append([])
    profiles_Bx.append([])
    profiles_By.append([])
    profiles_Bz.append([])
    profiles_Bl.append([])
    profiles_Bm.append([])
    profiles_Bn.append([])
    profiles_V.append([])
    profiles_Vx.append([])
    profiles_Vy.append([])
    profiles_Vz.append([])
    profiles_rho.append([])
    profiles_T.append([])
    vsc_eigenvalues.append([])
    vsc_eigenvectors.append([])

if len(sys.argv)==4:
    start= int(sys.argv[1])
    stop = int(sys.argv[2])
    cadence = int(sys.argv[3])
    figureTime = 0
elif len(sys.argv) == 5:
    start= int(sys.argv[1])
    stop = int(sys.argv[2])
    cadence = int(sys.argv[3])
    figureTime = int(sys.argv[4])
else:
    sys.stderr.write("Usage: plot_VSC.py <starting_index> <final_index+1> <cadence> \n")
    sys.stderr.write("OR\n")
    sys.stderr.write("Usage: plot_VSC.py <starting_index> <final_index+1> <cadence> <figure time>\n")
    sys.exit()


inputLocation="/wrk-vakka/group/spacephysics/vlasiator/3D/FHA/bulk1/"
outputLocation="./"
outputfile = outputLocation + "VSC_plot.png"

times=np.arange(start,stop,cadence)

for i,time in enumerate(times):
    
    timefull=str(time).rjust(7, '0')
    file_name=inputLocation+"bulk1."+timefull+".vlsv"

    print("Extracting timestep "+str(time)+" s")
    f=pt.vlsvfile.VlsvReader(file_name=file_name)

    profiles_times.append(f.read_parameter("time"))
    
    for sc in range(numVSC):
        cellid = f.get_cellid(pos[sc])

        B = f.read_variable("vg_b_vol", cellids=cellid)
        profiles_B[sc].append(np.linalg.norm(B)*1e9)
        profiles_Bx[sc].append(B[0]*1e9)
        profiles_By[sc].append(B[1]*1e9)
        profiles_Bz[sc].append(B[2]*1e9)


### MVA after Lucile Turc's code
def minimum_variance_analysis_Vlasiator(Bx,By,Bz):
    #______________________________ MINIMUM VARIANCE ANALYSIS ____________________________________

    # Creating two matrices of zeros, 'mag0' & 'mag'
    # number of rows = length of magnetic field list
    # number of columns = 3 (Bx, By, Bz)

    ndata = len(Bx)
    mag0 = np.zeros((ndata,3))
    mag = np.zeros((ndata,3))

    # mag0 is the matrix 3 column matrix of the x, y and z magnetic field values
    # mag is the mean adjusted magnetic field matrix
    mag0[:,0] = Bx
    mag0[:,1] = By
    mag0[:,2] = Bz

    mag[:,0] = np.copy(Bx) - ((1/float(ndata))*(sum(Bx)))
    mag[:,1] = np.copy(By) - ((1/float(ndata))*(sum(By)))
    mag[:,2] = np.copy(Bz) - ((1/float(ndata))*(sum(Bz)))

    # calculating the covariance matrix, Q

    magT = np.transpose(mag)
    Q = (np.dot(magT,mag))/ndata

    I = np.identity(3)

    # calculating the eigenvalues and eigenvectors
    eigenvalues, eigenvectors = np.linalg.eig(Q)

    # sorts the eigenvalues and eigenvector in descending size
    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:,idx]

    # redefines the last column of the eigenvector as the wave vector
    wave_vector = np.zeros((3))
    for i in range(0,3):
        wave_vector[i] = eigenvectors[i,2]

    # calculating the ratio of the eigenvectors
    ratio1 = eigenvalues[0]/eigenvalues[1]
    ratio2 = eigenvalues[0]/eigenvalues[2]
    ratio3 = eigenvalues[1]/eigenvalues[2]

    # tranforming the magnetic field by multiplying by the eigenvector matrix
    Bfield = np.array([(Bx[i],By[i],Bz[i]) for i in range(0,len(Bx))])
    Y = np.dot(Bfield,eigenvectors)
    eigenvectors = np.transpose(eigenvectors)

    # calculating the error
    error = math.sqrt((eigenvalues[1]*eigenvalues[2])/((ndata-1)*
        (eigenvalues[1]-eigenvalues[2])*(eigenvalues[1]-eigenvalues[2])))*(180/math.pi)

    # single value decomposition
    V, U, W = np.linalg.svd(mag)


    Bl = Y[:,0]
    Bm = Y[:,1]
    Bn = Y[:,2]

    print(f'Eigenvalues: {eigenvalues[0]:.2f} {eigenvalues[1]:.2f} {eigenvalues[2]:.2f}')
    print(f'{eigenvectors[0,0]:.2f} {eigenvectors[0,1]:.2f} {eigenvectors[0,2]:.2f}')
    print(f'{eigenvectors[1,0]:.2f} {eigenvectors[1,1]:.2f} {eigenvectors[1,2]:.2f}')
    print(f'{eigenvectors[2,0]:.2f} {eigenvectors[2,1]:.2f} {eigenvectors[2,2]:.2f}')
    print(f'lambda_1/lambda_2 = {ratio1:.2f}')
    print(f'lambda_1/lambda_3 = {ratio2:.2f}')
    print(f'lambda_2/lambda_3 = {ratio3:.2f}')

    return Bl, Bm, Bn, eigenvalues, eigenvectors


for sc in range(numVSC):
    profiles_Bl[sc], profiles_Bm[sc], profiles_Bn[sc], vsc_eigenvalues[sc], vsc_eigenvectors[sc] = minimum_variance_analysis_Vlasiator(profiles_Bx[sc], profiles_By[sc], profiles_Bz[sc])
    profiles_Blmnmag[sc] = [math.sqrt(profiles_Bl[sc][i]**2 + profiles_Bm[sc][i]**2 + profiles_Bn[sc][i]**2) for i in range(len(profiles_Bl[sc]))]

vsc_eigenvectors = np.array(vsc_eigenvectors)
vsc_eigenvalues = np.array(vsc_eigenvalues)

### end of MVA magic

pos=np.array(pos)

gsvunit=6
gshunitprofile=9
gsbuf=1
inter=0
cbwid=3
dpi=150
matplotlib.rcParams.update({'font.size': 35})

for i,time in enumerate(times):
    if figureTime != 0 and figureTime != time:
        continue
    timefull=str(time).rjust(7, '0')
    file_name=inputLocation+"bulk1."+timefull+".vlsv"
    f=pt.vlsvfile.VlsvReader(file_name=file_name)

    print("Slicing timestep "+str(time)+" s")
    # Init figure    
    fig = pl.figure(figsize=(2*gsvunit+gshunitprofile+cbwid+gsbuf+inter,4*gsvunit+3*gsbuf), dpi=dpi)
    gs = fig.add_gridspec(4*gsvunit+3*gsbuf, 2*gsvunit+gshunitprofile+cbwid+gsbuf+inter)

    fig.add_subplot(gs[:gsvunit,                            2*gsvunit+cbwid+gsbuf+inter:2*gsvunit+cbwid+gsbuf+inter+gshunitprofile]) # VSC 0
    fig.add_subplot(gs[gsvunit+gsbuf:2*gsvunit+gsbuf,       2*gsvunit+cbwid+gsbuf+inter:2*gsvunit+cbwid+gsbuf+inter+gshunitprofile]) # VSC 1
    fig.add_subplot(gs[2*gsvunit+2*gsbuf:3*gsvunit+2*gsbuf, 2*gsvunit+cbwid+gsbuf+inter:2*gsvunit+cbwid+gsbuf+inter+gshunitprofile]) # VSC 2
    fig.add_subplot(gs[3*gsvunit+3*gsbuf:,                  2*gsvunit+cbwid+gsbuf+inter:2*gsvunit+cbwid+gsbuf+inter+gshunitprofile]) # VSC 3
    fig.add_subplot(gs[:2*gsvunit, :2*gsvunit])                      # slice top
    fig.add_subplot(gs[2*gsvunit+3*gsbuf:, :2*gsvunit])                      # slice bottom
    fig.add_subplot(gs[:2*gsvunit, 2*gsvunit:2*gsvunit+cbwid//2]) # colorbar top
    fig.add_subplot(gs[2*gsvunit+3*gsbuf:, 2*gsvunit:2*gsvunit+cbwid//2]) # colorbar right

    axes = fig.get_axes()
    for ax in axes[0:3]:
#        for label in ax.xaxis.get_ticklabels():
#            label.set_visible(False)
        ax.set_xlabel("")
    axes[3].set_xlabel(r"Simulation time [s]")

    # plot each profile
    for sc in range(numVSC):
        ax=axes[sc]
        ax.set(ylim=(-13 , 13))
        ax.grid(linewidth=lw, color=colorseq[sc], alpha=0.2)
        for spine in ('top', 'bottom', 'right', 'left'):
            ax.spines[spine].set_color(colorseq[sc])
            ax.spines[spine].set_linewidth(3)
        ax.tick_params(width=3, length=9, color=colorseq[sc])
        #ax.set_ylabel(r"$B$ [nT]", fontsize=24)
        if plot_lmn == False:
            ax.plot(profiles_times,profiles_B[sc], lw=lw, color='k', ls='solid', label="$B$")
            ax.plot(profiles_times,profiles_Bx[sc], lw=lw, color=bgrseq[0], ls=linestyleseq[0], label="$B_x$")
            ax.plot(profiles_times,profiles_By[sc], lw=lw, color=bgrseq[1], ls=linestyleseq[1], label="$B_y$")
            ax.plot(profiles_times,profiles_Bz[sc], lw=lw, color=bgrseq[2], ls=linestyleseq[2], label="$B_z\,\mathrm{[nT]}$")
        else:
            # LMN from MVA
            ax.plot(profiles_times,profiles_Blmnmag[sc], lw=lw, color='k', ls='solid', label="$B$")
            ax.plot(profiles_times,profiles_Bl[sc], lw=lw, color=bgrseq[0], ls=linestyleseq[0], label="$B_l$")
            ax.plot(profiles_times,profiles_Bm[sc], lw=lw, color=bgrseq[1], ls=linestyleseq[1], label="$B_m$")
            ax.plot(profiles_times,profiles_Bn[sc], lw=lw, color=bgrseq[2], ls=linestyleseq[2], label="$B_n\,\mathrm{[nT]}$")
            # LMN from circle proxy
            #ax.plot(profiles_times,profiles_B[sc], lw=lw, color='k', ls='solid', label="$B$")
            #ax.plot(profiles_times,profiles_Bx[sc], lw=lw, color='blue', ls='dashed', label="$B_x$")
            #ax.plot(profiles_times,np.array(profiles_Bx[sc])*np.cos(angle[sc]) - np.array(profiles_By[sc])*np.sin(angle[sc]), lw=lw, color='green', ls='dashdot', label="$B_n$")
            #ax.plot(profiles_times,np.array(profiles_Bx[sc])*np.sin(angle[sc]) + np.array(profiles_By[sc])*np.cos(angle[sc]), lw=lw, color='red', ls='dotted', label="$B_l\,\mathrm{[nT]}$")
        ax.axvline(time, lw=lw, color='k')
        if figureTime != 0:
            panel_label_bbox = dict(boxstyle='square,pad=0.2', fc='white', ec='black', alpha=0.85, lw=thick)
            ax.text(stop-1.4, 11, panel_labels[sc], fontsize=42, bbox=panel_label_bbox, ha="right", va="top")
            for i in range(3):
                panel_label_bbox = dict(boxstyle='square,pad=0.2', fc='white', ec=bgrseq[i], alpha=0.85, lw=thick)
                ax.text(stop+5,12.5-i*8,"$\lambda_"+lmn[i]+"=$"+f"{vsc_eigenvalues[sc,i]:.2f}\n"+"$\mathbf{u}_"+lmn[i]+"=($"+f"{vsc_eigenvectors[sc,i,0]:.2f},{vsc_eigenvectors[sc,i,1]:.2f},{vsc_eigenvectors[sc,i,2]:.2f}"+"$)$", fontsize=32, color=bgrseq[i], bbox=panel_label_bbox, ha="left", va="top")

        if(sc==0):
            leg = ax.legend(loc='lower right', bbox_to_anchor=(1,1.08), ncol=4, frameon=True, edgecolor='black', handlelength=1.5, handletextpad=0.2, columnspacing=0.5, borderaxespad=0.1, borderpad=0.5, fancybox=False)
            leg.get_frame().set_linewidth(thick)

    # plot slice
    ax=axes[numVSC]
    if figureTime != 0:
        panel_label_bbox = dict(boxstyle='square,pad=0.2', fc='white', ec='black', alpha=0.85, lw=thick)
        ax.text(-20.5, -0.5, panel_labels[numVSC], fontsize=42, bbox=panel_label_bbox, ha="left", va="top")
    cbax=axes[numVSC+2]
    pt.plot.plot_colormap3dslice(
        circle=circleradius,circlex=circlex,circley=circley,
        vlsvobj=f,
        var="vg_b_vol",
#        external=[pt.plot.plot_helpers.reflevels, pt.plot.plot_helpers.reflevels_aligned],
#        fsaved="black",
#        fsavedlinewidth=3,
#        amr=[1,2,3],
#        amrlinewidths=[4,5,3],
#        amrcolours=["cyan", "magenta", "lime"],
#        amrlinestyles=["dashed", "dotted", "dashdot"],
        fluxrope=contourvalues,
        fluxropecolour=contourcolors,
        fluxropelinewidth=contourlinewidths,
        fluxropelinestyle="solid",
        operator="x", cbtitle="$B_x\,\mathrm{[nT]}$", boxre=[-21,-11,-10,0], normal='x', colormap='coolwarm',vmin=-100,vmax=100, symlog=True, usesci=None, vscale=1e9, step=time, streamlines = "vg_b_vol", streamlinedensity=2, streamlinethick=3, streamlinecolor='k', cutpointre=-10, axes=ax, cbaxes=cbax, scale=4, thick=thick)
#        operator="x", normal='x', colormap='bwr',vmin=-1e-7,vmax=1e-7, symlog=True, step=time, streamlines = "vg_b_vol", streamlinedensity=2, streamlinethick=3, streamlinecolor='k', cutpointre=-10, axes=ax, cbaxes=cbax, scale=3)
    ax.scatter(pos[:,1]/Re, pos[:,2]/Re, color=colorseq, s=1200, marker="o")
#    ax.quiver(pos[:,1]/Re, pos[:,2]/Re, np.ones(numVSC), np.ones(numVSC), angles=np.rad2deg(angle), color=colorseq, scale=10)
#    ax.quiver(pos[:,1]/Re, pos[:,2]/Re, np.ones(numVSC), np.ones(numVSC), angles=np.rad2deg(angle)-90, color=colorseq, scale=10)
    for i in range(3):
        for v in range(numVSC):
            headwidth = headwidth_default
            headlength = headlength_default
            headaxislength = headaxislength_default
            if vsc_eigenvectors[v,i,0] < 0:
                headwidth = 0
                headlength = 0
                headaxislength = 0
            ax.quiver(pos[v,1]/Re, pos[v,2]/Re, np.array(vsc_eigenvectors)[v,i,1], np.array(vsc_eigenvectors)[v,i,2], color=bgrseq[i], linestyle=linestyleseq[i], width=arrowwidth, scale=8, headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, zorder=triad_zorder)

    ax.text(-18.5, 1.3, f'            Into / out of plane', fontsize=36, bbox=dict(facecolor=(0,0,0,0), edgecolor='black', boxstyle='square,pad=0.2', lw=thick), ha='left', va='top')
    quiv = ax.quiver(100, 100, 101, 101, color='black', pivot='tip', scale=5, width=arrowwidth, headwidth=0, headlength=0, headaxislength=0)
    qk = ax.quiverkey(quiv, X=0.3, Y=1.105, U=0.5, angle=30, label=' ', labelpos='E', zorder=10)
    qk.pivot['E']='middle'
    quiv = ax.quiver(100, 100, 101, 101, color='black', pivot='tail', scale=5, width=arrowwidth, headwidth=headwidth_default, headlength=headlength_default, headaxislength=headaxislength_default)
    qk = ax.quiverkey(quiv, X=0.4, Y=1.105, U=0.5, angle=30, label=' ', labelpos='E', zorder=10)
    qk.pivot['E']='middle'



#    pt.plot.plot_colormap3dslice(vlsvobj=f, var="vg_j", operator="x", boxre=[-20,-10,-10,0], normal='x', colormap='bwr', vmin=-3e-9, vmax=3e-9, lin=9, usesci=True, step=time, streamlines = "vg_b_vol", streamlinedensity=2, streamlinethick=3, streamlinecolor='k', cutpointre=-10, axes=ax, cbaxes=cbax, cbtitle="fluxrope", scale=3)
    ax.xaxis.label.set_size(42)
    ax.yaxis.label.set_size(42)
    ax=axes[numVSC+1]
    if figureTime != 0:
        panel_label_bbox = dict(boxstyle='square,pad=0.2', fc='white', ec='black', alpha=0.85, lw=thick)
        ax.text(-14.5, -0.5, panel_labels[numVSC+1], fontsize=42, bbox=panel_label_bbox, ha="left", va="top")
    cbax=axes[numVSC+3]
    pt.plot.plot_colormap3dslice(
        vlsvobj=f,
        var="vg_b_vol",
#        external=[pt.plot.plot_helpers.reflevels, pt.plot.plot_helpers.reflevels_aligned],
#        fsaved="black",
#        fsavedlinewidth=3,
#        amr=[1,2,3],
#        amrlinewidths=[4,5,3],
#        amrcolours=["cyan", "magenta", "lime"],
#        amrlinestyles=["dashed", "dotted", "dashdot"],
        fluxrope=contourvalues,
        fluxropecolour=contourcolors,
        fluxropelinewidth=contourlinewidths,
        fluxropelinestyle="solid",
        operator="y", cbtitle="$B_y\,\mathrm{[nT]}", boxre=[-15,-5,-10,0], normal='y', colormap='coolwarm',vmin=-100,vmax=100, symlog=True, usesci=None, vscale=1e9, step=time, streamlines = "vg_b_vol", streamlinedensity=2, streamlinethick=3, streamlinecolor='k', cutpointre=-17,  axes=ax, cbaxes=cbax, scale=4, thick=thick)
#        operator="y", normal='y', colormap='bwr',vmin=-1e-7,vmax=1e-7, symlog=True, step=time, streamlines = "vg_b_vol", streamlinedensity=2, streamlinethick=3, streamlinecolor='k', cutpointre=-17, axes=ax, cbaxes=cbax, scale=3)
    ax.scatter(pos[:,0]/Re, pos[:,2]/Re, color=colorseq, s=1200, marker="o")
    for i in range(3):
        for v in range(numVSC):
            headwidth = headwidth_default
            headlength = headlength_default
            headaxislength = headaxislength_default
            if vsc_eigenvectors[v,i,1] < 0:
                headwidth = 0
                headlength = 0
                headaxislength = 0
            ax.quiver(pos[v,0]/Re, pos[v,2]/Re, np.array(vsc_eigenvectors)[v,i,0], np.array(vsc_eigenvectors)[v,i,2], color=bgrseq[i], width=arrowwidth, scale=8, linestyle=linestyleseq[i], headwidth=headwidth, headlength=headlength, headaxislength=headaxislength, zorder=triad_zorder)

    ax.xaxis.label.set_size(42)
    ax.yaxis.label.set_size(42)

    time_label_bbox = dict(boxstyle='square,pad=0.1', fc='white', ec='black', alpha=0.85, lw=thick)
    fig.suptitle("t="+str(time)+" s", x=0.07, y=0.93, ha="left", va="top", fontsize=50, bbox=time_label_bbox)


    if figureTime != 0:
        if plot_lmn == False:
            fig.savefig("FHA_VSC_shrimp_slices_figure_"+timefull+".png", bbox_inches='tight')
        else:
            fig.savefig("FHA_VSC_shrimp_slices_LMN_figure_"+timefull+".png", bbox_inches='tight')
    else:
        if plot_lmn == False:
            fig.savefig("FHA_VSC_shrimp_slices_"+timefull+".png", bbox_inches='tight')
        else:
            fig.savefig("FHA_VSC_shrimp_slices_LMN_"+timefull+".png", bbox_inches='tight')
#    fig.savefig("FHA_VSC_shrimp_slice_jx_"+timefull+".png", bbox_inches='tight')

    pl.close()
