import analysator as pt
import sys, os
import numpy as np
import traceback
import inspect
import logging
import re
import argparse
from testpackage_template_maker import call_replace

def source_file_name(filename,fileLocation,time):
    if filename is None:
        if '2D' not in fileLocation:
            bulkname = "bulk1."+str(time).rjust(7,'0')+".vlsv"

        else:
            bulkname = "bulk."+str(time).rjust(7,'0')+".vlsv"
    else:
        bulkname=filename
    return bulkname


    
argp=argparse.ArgumentParser(
    prog='Analysator Testpackage',
    description='Outputs test plots'
)
argp.add_argument("jobcount",type=int)
argp.add_argument("jobindex",type=int)

argp.add_argument('funcs',type=str,help="function/list of functions to test, if none give does all.",nargs='*')
cmd_args=argp.parse_args()
funcs_to_use=cmd_args.funcs
if "pass" in funcs_to_use:
    quit()

datalocation = "/wrk-vakka/group/spacephysics/vlasiator"
runs = []

#Change this to make it produce same plots as testpackage_vdf and testpackage_colormap used to do

runs.append( { 'name': 'ABC',
                 'verifydir': '/ABC/', 
                 'fileLocation': datalocation+'/2D/ABC/bulk/',
                 'fluxLocation': datalocation+'/2D/ABC/flux/',
                'skipped_args':{'plot_vdf':{'step':''},'plot_vdf_profiles':{'bpara':'','bperp':'','step':''},
                                'plot_vdfdiff':{'filedir':'','bpara':'','bperp':''}},
                 'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles','plot_vdfdiff'],
                 'pops': ['avgs'],
                 'time': 1000,
                 'singletime': False,
                 'filename': None,
                 'nosubpops': False, # backstreaming / non-backstreaming
                 'vlasiator5': False,
                 'cavitonparams': [6.6e6,2.64e6,4.e-9,10] } )

runs.append( { 'name': 'BED', 
                 'verifydir': '/BED/', 
                'fluxLocation': None,
                 'fileLocation': datalocation+'/2D/BED/bulk/',
                 'pops': ['avgs'],
                 'time': 2000,
                 'singletime':True,
                'nosubpops': False, # backstreaming / non-backstreaming
                 'vlasiator5': False,
                 'skipped_args':{'plot_vdf_profiles':{'bpara':'','bperp':'','step':''}},
                 'funcs': ['plot_vdf','plot_vdf_profiles'],
                 'filename': None ,
                'cavitonparams': [6.6e6,2.64e6,4.e-9,10]
                })

runs.append( { 'name': 'FID',
                 'verifydir': '/FID/', 
                 'fileLocation': datalocation+'/3D/FID/bulk1/',
                 'fluxLocation': None,
                 'funcs': ['plot_colormap3dslice','plot_ionosphere','plot_isosurface'],
                 'pops': ['avgs'],
                 'time': 1000,
                'skipped_args':{'plot_ionosphere':{"var":["ig_z","ig_p","ig_source","ig_residual"]},'ALL':{'expression':''}}, #ig_zz and ig_pp is also skipped on purpose
                 'singletime': False,
                 'filename': None, #restart file
                 'manualcall':False,
                 'nosubpops': True, # backstreaming / non-backstreaming
                 'vlasiator5': True,
                 'cavitonparams': [6.6e6,2.64e6,4.e-9,10] } )

runs.append( { 'name': 'FHA',
                 'verifydir': '/FHA/', 
                 'fileLocation': datalocation+'/3D/FHA/bulk1/',
                 'fluxLocation': None,
                 'funcs': ['plot_colormap3dslice','plot_ionosphere','plot_isosurface','plot_vdf','plot_vdfdiff'],
                 'pops': ['avgs'],
                 'time': 1000,
                'skipped_args':{'ALL':{'expression':''}},
                 'singletime': False,
                 'filename': None, #restart file
                 'manualcall':False,
                 'nosubpops': True, # backstreaming / non-backstreaming
                 'vlasiator5': True,
                 'cavitonparams': [6.6e6,2.64e6,4.e-9,10] } )

runs.append( { 'name': 'BCQ',
                'verifydir': '/BCQ/', 
                'fileLocation': datalocation+'/2D/BCQ/bulk/',
                'fluxLocation': None,
                'singletime': False,
                'pops': ['avgs'],
                'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles','plot_vdfdiff'],
                'time': 2000,
                 'skipped_args':{'plot_vdf_profiles':{'bpara':'','bperp':'','step':''}},
                'filename': None,
                'vlasiator5':False,
                'nosubpops':False,
                'cavitonparams': [2.0e6,0.8e6,4.e-9,10] 
                  } )

runs.append( { 'name': 'BCQr',
                 'verifydir': '/BCQr/', 
                 'funcs': ['plot_vdf','plot_vdf_profiles'],
                 'fileLocation': datalocation+'/2D/BCQ/restart/',
                 'pops': ['avgs'],
                'fluxLocation': None,
                 'singletime': True, # neighboring bulk files not available
                 'time': 1361,
                'skipped_args':{'plot_vdf':{'step':''}},
                 'manualcall':False,
                 'vlasiator5': False,
                 'nosubpops': False, # thermal / non-thermal
                 'filename': 'restart.0001361.vlsv',
                'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )

runs.append( { 'name': 'BGA',
                 'verifydir': '/BGA/', 
                 'fileLocation': datalocation+'/2D/BGA/zero_ehall_layers_23/',
                 'fluxLocation': None,
                 'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles','plot_vdfdiff'],
                 'pops': ['proton'],
                 'skipped_args':None,
                 'time': 380,
                 'manualcall':False,
                 'singletime': True, # neighboring bulk files not available
                 'filename': None,
                 'vlasiator5': True,
                 'nosubpops': True, # thermal / non-thermal
                 'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )
                     


runs.append( { 'name': 'BFDr',
                 'verifydir': '/BFDr/', 
                'funcs': ['plot_vdf'],
                 'fileLocation': datalocation+'/2D/BFD/restart/',
                 'pops': ['avgs'],
                'fluxLocation': None,
                 'singletime': True, # neighboring bulk files not available
                 'time': 1126,
                'manualcall':False,
                'skipped_args':{'plot_vdf':{'step':''}},
                'nosubpops': False, # thermal / non-thermal
                'vlasiator5': False,
                 'filename': 'restart.0001126.vlsv',
                'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )

runs.append( { 'name': 'BFD',
                 'verifydir': '/BFD/', 
                 'fileLocation': datalocation+'/2D/BFD/bulk/',
                 'fluxLocation': datalocation+'/2D/BFD/fluxfunction/',
                 'fluxprefix': 'bulk.',
                 'funcs': ['plot_colormap','plot_vdf','plot_vdf_profiles'],
                 'skipped_args':None,
                 'pops': ['proton','helium'],
                 'time': 2000,
                 'singletime': False,
                 'filename': None,
                  'manualcall':False,
                 'nosubpops': False, # backstreaming / non-backstreaming
                 'vlasiator5': False,
                 'cavitonparams': [2.0e6,0.8e6,4.e-9,10] } )

#First arg in tuple is the one that needs to be there, second etc can be there without the first
required_args ={
    "plot_vdf":[(["coordre","coordinates","cellids"],["coordre=REPLACECOORDRE"]),([("filedir","step"),'vlsvobj','filename'],None)],
    "plot_vdf_profiles":[(["coordre","coordinates","cellids"],["coordre=REPLACECOORDRE"]),([("filedir","step"),'vlsvobj','filename'],None)],
    "plot_vdfdiff":[(["coordre","coordinates","cellids"],["coordre=REPLACECOORDRE"]),([("filedir","step"),'vlsvobj','filename'],None)],
    "plot_isosurface":[([("surf_step","surf_var")],["surf_step=10","surf_var='vg_rho'"]),([("filedir","step"),'vlsvobj','filename'],None)]
    
}

calls = []
callrunids = []
callrunindex = []
funcids=[]

for i,run in enumerate(runs):
    # bulk and restart files
    vlasiator5 = run['vlasiator5']
    filename = run['filename']
    fileLocation = run['fileLocation']
    singletime = run['singletime']
    nosubpops = run['nosubpops']
    fluxLocation = run['fluxLocation']

#    if run['name']!='FID':
#        continue

    if not funcs_to_use:
        functions = run['funcs']

    else:
        functions = sorted(list(set(run['funcs']) & set(funcs_to_use)))
        if not functions:
            continue


    for j,func in enumerate(functions):
        callindex = 0
        
        #try to import the list of calls corresponding to the function to be tested. Skip if not found
        try:
            exec(f'import testpackage_definitions.testpackage_{func} as testpackage_{func}')
        except:
            raise IOError(f"testpackage_{func} could not be imported, check that the file exists and is in the correct folder")


        #Get the list of calls from the imported file, set list to empty list if list not foud in the file
        for call_list in ["restartcalls","nonrestartcalls","multipopcalls","v5restartcalls","v5nonrestartcalls","v5multipopcalls"]:
            try:
                exec(f'{call_list}=testpackage_{func}.{call_list}')
            except:
                exec(f'{call_list}=[]')

        skipped_args=run['skipped_args']


        if filename is not None:
            calls_in=v5restartcalls if vlasiator5 else restartcalls
            for call in calls_in:
                if vlasiator5:
                    call = call.replace("var='vg_v'","var='vg_restart_v'")
                else:
                    call = call.replace("var='V'","var='restart_V'")
                if skipped_args:
                    call=call_replace(call,func,skipped_args,required_args)
                if call is not None: 
                    callrunids.append(i)
                    calls.append(call)
                    callrunindex.append(callindex)
                    callindex += 1
                    funcids.append(j)


        # non-restart files
        if filename is None:
            #non restart calls
            calls_in=v5nonrestartcalls if vlasiator5 else nonrestartcalls
            for call in calls_in:

                # Skip flux function calls if no flux files
                if ("fluxdir" in call or "fluxfile" in call) and fluxLocation is None:
                    continue

                #change the extrnal and expression calls to correct format
                if f"expression=testpackage_{func}" not in call and "expression=" in call:
                    call=call.replace("expression=",f"expression=testpackage_{func}.cexp.")
                if f"external=testpackage_{func}" not in call and "external=" in call:
                    call=call.replace("external=",f"external=testpackage_{func}.cexp.")

                # skip time integration if only one file available
                if "pass_times" in call and singletime:
                    continue
                # thermal / non-thermal subpopulations
                if vlasiator5 and (("_thermal" in call) or ("_nonthermal" in call)) and nosubpops:
                    continue
                elif (("_backstream" in call) or ("_nonbackstream" in call)) and nosubpops:
                    continue
                if skipped_args:
                    call=call_replace(call,func,skipped_args,required_args)
                if call is not None:
                    callrunids.append(i)
                    calls.append(call)
                    callrunindex.append(callindex)
                    callindex += 1
                    funcids.append(j)


                
            #multipop calls
            calls_in=v5multipopcalls if vlasiator5 else multipopcalls
            for pop in run['pops']:
                if pop != 'avgs':
                    for call in calls_in:
                        # Skip flux function calls if no flux files
                        if "flux" in call and fluxLocation is None:
                            continue
                        # skip time integration if only one file available
                        if "pass_times" in call and singletime:
                            continue
                        # thermal / non-thermal subpopulations
                        if vlasiator5 and (("_thermal" in call) or ("_nonthermal" in call)) and nosubpops:
                            continue
                        elif (("_backstream" in call) or ("_nonbackstream" in call)) and nosubpops:
                            continue
                        call = call.replace('REPLACEPOP',pop)
                        if skipped_args:
                            call=call_replace(call,func,skipped_args,required_args)
                        if call is not None:
                            callrunids.append(i)
                            calls.append(call)
                            callrunindex.append(callindex)
                            callindex += 1
                            funcids.append(j)



nteststot = len(callrunids)

# How many jobs? 
jobcount=cmd_args.jobcount
jobcurr=cmd_args.jobindex

increment = int(nteststot/jobcount)
remainder = nteststot - jobcount * increment
start=jobcurr * increment
end=start + increment
# Remainder frames are divvied out evenly among tasks
if jobcurr < remainder:
    start = start + jobcurr
    end = end + jobcurr + 1
else:
    start = start + remainder
    end = end + remainder



# Perform calls
for j in range(start,end):
    # Calculate which run
    jrun = callrunindex[j]
    runid = callrunids[j]
    runname = runs[runid]['name']

    call = calls[j]

    funcid=funcids[j] 
    

    if not funcs_to_use:
        func = runs[runid]['funcs'][funcid]

    else:        
        func = sorted(list(set(runs[runid]['funcs']) & set(funcs_to_use)))[funcid]
        if not func:
            continue

    verifydir = func+runs[runid]['verifydir']

    fileLocation = runs[runid]['fileLocation']
    fluxLocation = runs[runid]['fluxLocation']
    pops = runs[runid]['pops']
    time = runs[runid]['time']
    filename = runs[runid]['filename']
    vlasiator5 = runs[runid]['vlasiator5']
    singletime = runs[runid]['singletime']
    
    #set custom expression variables for plot_colormap
    if 'plot_colormap' in func:
        exec(f"testpackage_{func}.cexp.level_bow_shock = runs[runid]['cavitonparams'][0]")
        exec(f"testpackage_{func}.cexp.level_n_caviton = runs[runid]['cavitonparams'][1]")
        exec(f"testpackage_{func}.cexp.level_B_caviton = runs[runid]['cavitonparams'][2]")
        exec(f"testpackage_{func}.cexp.level_beta_SHFA = runs[runid]['cavitonparams'][3]")
    
        #Custom expression vlasiator5 toggle
        if vlasiator5:
            exec(f'testpackage_{func}.cexp.vlasiator5=True')
        else:
            exec(f'testpackage_{func}.cexp.vlasiator5=False')

    
    outputLocation=os.path.join(pt.plot.defaultoutputdir,verifydir)


    # Source data files
    bulkname = source_file_name(filename,fileLocation,time)

    if 'fluxprefix' in runs[runid]:
        fluxname = runs[runid]['fluxprefix']+str(time).rjust(7,'0')+".bin"
    else:
        fluxname = "flux."+str(time).rjust(7,'0')+".bin"

    if runs[runid]['name']=="ABC" and 'plot_vdf' in func:
        fileLocation=fileLocation.replace('bulk','distributions')
        bulkname = "distributions."+str(100).rjust(7,'0')+".vlsv" 
        time=100

    if func == 'plot_vdfdiff':
        time_offset=10
        bulkname_vdfdiff = source_file_name(filename,fileLocation,time-time_offset)
        call = call.replace('REPLACETIME2',"'"+str(time-time_offset)+"'")
        if runs[runid]['name']=="ABC":
            bulkname_vdfdiff=bulkname_vdfdiff.replace('bulk','distributions')
        g = pt.vlsvfile.VlsvReader(fileLocation+bulkname_vdfdiff)
    

    call = call.replace('REPLACEPREVINDEX',"'"+str(jrun-1).rjust(4,'0')+"'")
    call = call.replace('REPLACEINDEX',"'"+str(jrun).rjust(4,'0')+"'")
    call = call.replace('REPLACETIME',"'"+str(time)+"'")

    call = call.replace('REPLACECELLID','1')
    call = call.replace('REPLACECOORDRE','[10,0,0]')
    call = call.replace('REPLACECOORDINATES','[6.371e7,0,0]')
    call = call.replace('REPLACEMULTIPLECELLID','[1,51,101]')
    call = call.replace('REPLACEMULTIPLECOORDRE','[[10,0,0],[15,0,0],[20,0,0]]')
    call = call.replace('REPLACEMULTIPLECOORDINATES','[[6.371e7,0,0],[9.5565e7,0,0],[12.742e7,0,0]]')


    # Many different plots
    print(j, runid, jrun, call,fileLocation+bulkname)
    f = pt.vlsvfile.VlsvReader(fileLocation+bulkname)

    try:
        exec(call)

    except Exception as e:
        print("----------------------------\nFAILURE DURING CALL ",j, runname," \n```\n"+call+"```\n", repr(e))
        
        traceback.print_exc()
        print("END TRACE for call",j,"\n----------------------------")

