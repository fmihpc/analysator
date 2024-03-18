import sys,os,shutil,time
import numpy as np
import math

indx    = int(sys.argv[1])

pathRun = '/scratch/project_2000203/dubartma/dmumu/'

pathDmumuTemplate = "/users/dubartma/analysis/subgrid/TimeDependence/"
DmumuTemplate     = 'DmumuFit.sh'

with open(pathDmumuTemplate+DmumuTemplate,'r') as file:
    fileDmumu = file.read()

Taniso = np.linspace(np.log10(1.0),np.log10(10.0),20)
betaPara = np.linspace(np.log10(0.01),np.log10(100),20)

bulkStart  = int(sys.argv[2])
bulkEnd    = int(sys.argv[3])
nbins      = 30
interval   = int(sys.argv[4])
twindow    = int(sys.argv[5]) #Time average interval for f(mu)

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

       # Create solver file ---------------
       printDmumu = fileDmumu
       printDmumu = printDmumu.replace('$jobName$'   ,'Dmumu_'+newname)
       printDmumu = printDmumu.replace('$start$',str(bulkStart))
       printDmumu = printDmumu.replace('$end$'  ,str(bulkEnd))
       printDmumu = printDmumu.replace('$folder$'    ,newname)
       printDmumu = printDmumu.replace('$interval$'  ,str(interval))
       printDmumu = printDmumu.replace('$twindow$'   ,str(twindow))

       with open(pathDmumuTemplate + '/batch/Dmumu_'+newname+'.sh' ,'w') as file:
           file.write(printDmumu)

       os.chdir(pathDmumuTemplate + '/batch/')
       os.system("sbatch "+'Dmumu_'+newname+'.sh')
       os.chdir(pathDmumuTemplate)
