#Script which generates cfg and sh files and run local runs for diffusion
#Uses LossCone.cfg and LossCone.sh as templates
#Hardcoded with personal paths, so need to be changed
#Creates the needed path to be accessed by the scripts which compute Dmumu

import sys,os,shutil,time
import numpy as np
import math

indx = int(sys.argv[1]) # From 0 to 3, batch of 100 runs

pathRun = '/scratch/project_2000203/dubartma/dmumu/'

cfgTemplate = "LossCone.cfg"
shTemplate  = "LossCone.sh"

pathRunTemplate = '/scratch/project_2000203/dubartma/dmumu/'

with open(pathRunTemplate+cfgTemplate,'r') as file:
   filecfg = file.read()

with open(pathRunTemplate+shTemplate,'r') as file:
   filesh = file.read()

Taniso   = np.linspace(np.log10(1.0),np.log10(10.0),20)
betaPara = np.linspace(np.log10(0.01),np.log10(100),20)

n   = 3.0e6 # Density, similar to BCQ mid sheath
kB  = 1.380649e-23
mu0 = 1.25663706212e-6
q   = 1.602176634e-19 
m   = 1.67262192e-27

Tx = 5.72e6 # Parallel temperature (fixed) similar to BCQ mid sheath

Bx  = np.sqrt(n*kB*Tx*2.0*mu0/10**betaPara) # By,z = 0
fcp = q*Bx/(2.0*np.pi*m) # Proton cyclotron frequency

TimeOfGyration = 1.0/fcp # Time for proton to do one gyration

TimeMax  = 20*TimeOfGyration  # Max time of simulation to model 20 gyrations
dtOutput = TimeOfGyration/10  # Output data every 10th of a gyration
dtS      = TimeOfGyration/100 # Max Simulation time step

res       = 80000 # 80 km cells wave length is around 800km
cellNUM   = 100
lengthMAX = res*cellNUM # fits around 10 wavelengths in the box

for i in range(0,len(Taniso)):
    for j in range(indx*5,indx*5+5):

       newname = "300_T"+str(round(10**Taniso[i],3))+"_Beta"+str(round(10**betaPara[j],3))
       newpath = pathRun+newname
       pathData = "/scratch/project_2000203/dubartma/data/dmumu/"+newname+"/"
       pathFig  = "/scratch/project_2000203/dubartma/fig/dmumu/" +newname+"/"
       try:
           os.mkdir(pathData)
       except:
           pass
       try:
           os.mkdir(pathData+"/mu/")
       except:
           pass
       try:
           os.mkdir(pathData+"/Dmumu/")
       except:
          pass
       try:
           os.mkdir(pathData+"/variables/")
       except:
           pass
       try:
          os.mkdir(newpath)
       except:
          pass
       try:
          os.mkdir(newpath+"/bulk/")
       except:
          pass

       Tyz = Tx * 10**Taniso[i] # Perp temperature

       cfgName = newname+'.cfg'

       printcfg = filecfg
       printcfg = printcfg.replace('$dtOutput$'    ,"{:.3}".format(dtOutput[j]))
       printcfg = printcfg.replace('$cells$'       ,str(cellNUM))
       printcfg = printcfg.replace('$length$'      ,str(lengthMAX))
       printcfg = printcfg.replace('$res$'         ,str(res))
       printcfg = printcfg.replace('$tMax$'        ,"{:.3}".format(TimeMax[j]))
       printcfg = printcfg.replace('$dtS$'         ,"{:.3}".format(dtS[j]))
       printcfg = printcfg.replace('$BX0$'         ,"{:.2e}".format(Bx[j]))
       printcfg = printcfg.replace('$rho$'         ,"{:.2e}".format(n))
       printcfg = printcfg.replace('$TemperatureX$',"{:.2e}".format(Tx))
       printcfg = printcfg.replace('$TemperatureY$',"{:.2e}".format(Tyz))
       printcfg = printcfg.replace('$TemperatureZ$',"{:.2e}".format(Tyz))

       with open(newpath+"/"+cfgName, 'w') as file:
          file.write(printcfg)

       # Create Run sh file ---------------
       shName = newname+'.sh'

       printsh = filesh
       printsh = printsh.replace('$jobName$',newname)
       printsh = printsh.replace('$configfile$','"'+newpath+'/'+cfgName+'"')

       with open(newpath+"/"+shName, 'w') as file:
          file.write(printsh)

       os.chdir(newpath)
       os.system("sbatch "+shName)
       os.chdir(pathRunTemplate)

time.sleep(5)

