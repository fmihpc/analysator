import numpy as np
import sys,os
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
import matplotlib
import matplotlib.colors as colors
import scipy
import scipy.interpolate
import math
from scipy.interpolate import RegularGridInterpolator
import traceback

# Arguments
run       = sys.argv[1]
bulkStart = int(sys.argv[2])
bulkEnd   = int(sys.argv[3])
interval  = int(sys.argv[4])

rerun     = sys.argv[5]

# Paths
pathFig  = '/users/dubartma/analysis/fig/3D/'
pathLoad = '/users/dubartma/analysis/data/'

timetot = range(bulkStart, bulkEnd-2*interval+1,1)

betaPara = np.linspace(np.log10(0.01),np.log10(100),20)
Taniso   = np.linspace(np.log10(1.0),np.log10(10.0),20)

nu0Array    = np.empty(len(betaPara)*len(Taniso)*len(timetot))
betaArray   = np.empty(len(betaPara)*len(Taniso)*len(timetot))
TanisoArray = np.empty(len(betaPara)*len(Taniso)*len(timetot))
timeArray   = np.empty(len(betaPara)*len(Taniso)*len(timetot))

if rerun == 'True':

    print('Rerun')

    nu03D        = np.zeros([len(Taniso),len(betaPara),len(timetot)])

    for i in range(0,len(Taniso)):
        for j in range(0,len(betaPara)):
    
            newname   = "300_T"+str(round(10**Taniso[i],3))+"_Beta"+str(round(10**betaPara[j],3))
            pathData  = '/scratch/project_2000203/dubartma/data/'+run+'/'+newname+'/'

            try:

                nu03D[j,i,:] = np.load(pathData+'Dmumu/nu0Mean_'+run+'_'+newname+'_'+str(interval)+'.npy')            
                #nu03D[j,i,:] = np.load(pathData+'Dmumu/nu0Analytic_'+run+'_'+newname+'_'+str(interval)+'.npy')            

            except Exception as e:
                print(e)
                continue

    np.save(pathLoad+'nu03D.npy',nu03D)

else:
    print('Loading data')
    nu03D = np.load(pathLoad+'nu03D.npy')

TanisoAIC = 1 + 0.43/(10**betaPara + 0.0004)**(0.42)
TanisoMM  = 1 + 0.77/(10**betaPara + 0.016)**(0.76)

for t in timetot:

    BetaPara_curve = np.arange(0.01,100.01,0.01)
    Taniso_AIC     = 1 + 0.43/(BetaPara_curve + 0.0004)**(0.42) # PCI marginal threshold condition, Hellinger et al. 2006, Yoon et al. 2012
    Taniso_MM      = 1 + 0.77/(BetaPara_curve + 0.016)**(0.76)  # MM marginal threshold condition, Hellinger et al. 2006, Yoon et al. 2012
    
    # ----------- nu03D --------------------------------
  
    nu03D[:,0,:] = 0.0

    nu0plot = nu03D
    nu0min  = 5e-3

    for i in range(0,len(Taniso)):
        for j in range(0,len(betaPara)):
   
            if nu03D[j,i,t] < nu0min: nu0plot[j,i,t] = nu0min

            nu0Array   [t*len(betaPara)*len(Taniso)+i*len(betaPara)+j] = nu03D[j,i,t]
            TanisoArray[t*len(betaPara)*len(Taniso)+i*len(betaPara)+j] = 10**Taniso[i]
            betaArray  [t*len(betaPara)*len(Taniso)+i*len(betaPara)+j] = 10**betaPara[j]
            timeArray  [t*len(betaPara)*len(Taniso)+i*len(betaPara)+j] = (t+1)*0.1
  
    plt.pcolormesh(10**betaPara,10**Taniso,nu0plot[:,:,t].T,shading='nearest',norm=colors.LogNorm(vmin=nu0min,vmax=0.5),cmap='plasma')
    #plt.pcolormesh(10**betaPara,10**Taniso,nu03D[:,:,t].T,shading='nearest',vmin=5.0e-2,vmax=0.5,cmap='plasma')
    #plt.scatter(betaArray[t*len(betaPara)*len(Taniso):(t+1)*len(betaPara)*len(Taniso)],TanisoArray[t*len(betaPara)*len(Taniso):(t+1)*len(betaPara)*len(Taniso)],c=nu0Array[t*len(betaPara)*len(Taniso):(t+1)*len(betaPara)*len(Taniso)],vmin=0.05,vmax=0.5,cmap='plasma')

    plt.loglog(BetaPara_curve,Taniso_AIC,linewidth=2,color='white',label='PCI')
    plt.loglog(BetaPara_curve,Taniso_MM ,linewidth=2,color='white',linestyle='--',label='Mirror')
    
    plt.xlim(BetaPara_curve[0],BetaPara_curve[-1])
    plt.ylim(1.0,10.0)
    
    plt.xlabel('$\\beta_{\\parallel,0}$'           ,fontsize=15)
    plt.ylabel('${T_\\bot / T_\\parallel}_0$',fontsize=15)
    
    ax = plt.gca()
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    ax.yaxis.set_minor_formatter(FormatStrFormatter('%.0f'))
    
    plt.legend(loc='lower left')
    
    plt.colorbar()
    plt.text(200,11,'$\\nu_0$',fontsize=20,weight='bold')
    
    plt.suptitle('t = '+str(round((t+1.5)*0.1,1))+' $f_\\mathrm{cp}^{-1}$',fontsize=15)
    plt.savefig(pathFig+'Dmumu3D_'+str(t).rjust(7,'0')+'.png',dpi=300)
    print      (pathFig+'Dmumu3D_'+str(t).rjust(7,'0')+'.png')
    plt.close()
    
header = 't*fcp betaPara Taniso nu0'
data   = np.column_stack([timeArray,betaArray,TanisoArray,nu0Array])
np.savetxt(pathLoad+'nu03Ddt.dat',data,header=header)
print(pathLoad+'nu03Ddt.dat')

