import sys,os,shutil,time
import numpy as np
import math

indx    = int(sys.argv[1])

pathRun = '/scratch/project_2000203/dubartma/dmumu/'

pathDmumuTemplate = "/users/dubartma/analysis/subgrid/TimeDependence/"
VariablesTemplate = 'variables.sh'

with open(pathDmumuTemplate + VariablesTemplate,'r') as file:
    fileVariables = file.read()

Taniso = np.linspace(np.log10(1.0),np.log10(10.0),20)

betaPara = np.linspace(np.log10(0.01),np.log10(100),20)

bulkStart  = int(sys.argv[2])
bulkMax    = int(sys.argv[3])
nbins      = 30               #bins to build 1D mu space
interval   = int(sys.argv[4]) #interval to compute Dmumu (compute de difference in VDF between f(i) and f(i+interval). Has to be >=1

for i in range(0,len(Taniso)):
    for j in range(indx*5,indx*5+5):


       newname = "300_T"+str(round(10**Taniso[i],3))+"_Beta"+str(round(10**betaPara[j],3))
       newpath = pathRun+newname
       pathData = "/scratch/project_2000203/dubartma/data/dmumu/"+newname+"/"
       try:
           os.mkdir(pathDmumuTemplate+'batch/')
       except:
          pass
       try:
           os.mkdir(pathData)
       except:
          pass

       # Create variables sh file --------------
       printVariables = fileVariables
       printVariables = printVariables.replace('$jobName$'   ,'Variables_'+newname)
       printVariables = printVariables.replace('$start$'  ,str(bulkStart))
       printVariables = printVariables.replace('$end$'    ,str(bulkMax))
       printVariables = printVariables.replace('$folder$'    ,newname)
       printVariables = printVariables.replace('$interval$'  ,str(interval))

       with open(pathDmumuTemplate+"/batch/Variables_"+newname+'.sh', 'w') as file:
          file.write(printVariables)

       os.chdir(pathDmumuTemplate + '/batch/')
       os.system("sbatch "+"Variables_"+newname+'.sh')
       os.chdir(pathDmumuTemplate)
