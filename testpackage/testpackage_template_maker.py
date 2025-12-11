import analysator as pt
import sys, os
import numpy as np
import traceback
import inspect
import logging
import re
import argparse

required_args=False
def call_replace(call,func,skipped_args,required_args=required_args):
    #This is kind of scuffed maybe

    call=call.replace('REPLACEFUNC',func)

    #Get the arguments of the call
    args = re.search(r'\((.+)\)',call).group(1)
    args = [x[0] or x[1] or x[2] for x in re.findall(r'(\w+=[^,(\[]+)|(\w+=\(.+\))|(\w+=\[.+?\])',args)]
    named_parameters=[arg.split("=")[0] for arg in args]
    #Get parameters of the func
    function_pars=inspect.getfullargspec(eval("pt.plot."+func)).args
    #Remove args that are not present as parameters for the func
    args_out=[]

    #check that all required func args are set
    if required_args and func in required_args.keys():
        for required_tuple in required_args[func]:
            required_params=required_tuple[0]
            default_params=required_tuple[1]

            check=False
            for param in required_params:
                if type(param)==tuple and param[0] not in named_parameters:
                    check=True
                elif any((all(r in named_parameters for r in param),(param in named_parameters))):
                    if not skipped_args or func not in skipped_args or not any(r in skipped_args[func].keys() for r in param):
                        check=True
                        break
            if not check:
                #Add parameters if there are default_params
                if default_params:
                    if type(default_params)==str:
                        default_params=[default_params]
                    for param in default_params:
                        if param not in named_parameters:
                            args_out.append(param)
                else:
                    return None


    #skip args if there are skipped args and append if called arg in function_pars
    for arg in args:
        if arg:
            if skipped_args:
                if func in skipped_args.keys():
                    skipped_args_dict=skipped_args[func]
                elif 'ALL' in skipped_args.keys():
                    skipped_args_dict=skipped_args['ALL']
                else:
                    skipped_args_dict=False
                call_args=arg.split("=")
                if skipped_args_dict and call_args[0] in skipped_args_dict.keys():
                    if type(skipped_args_dict[call_args[0]])==str and skipped_args_dict[call_args[0]] in call_args[1]:
                        continue
                    elif type(skipped_args_dict[call_args[0]])==list and any(arg_skip in call_args[1] for arg_skip in skipped_args_dict[call_args[0]]):  #list of args in dict value means OR  ex. {'var':["vg_rho","vg_phi"]}
                        continue
            if arg.split("=")[0] in function_pars:
                args_out.append(arg)
            #else:
            #    logging.warning(f"Argument {arg} removed from call {call}")
    
    if not args_out:
        return None
    args_out=filter(None,args_out)
    call=call[:call.rfind("(")+1]+",".join(args_out)+")"

    return call


if __name__=='__main__':


    argp=argparse.ArgumentParser(
        prog='Analysator Testpackage',
        description='Outputs test plots'
    )


    argp.add_argument('funcs',type=str,help="function/list of functions to test, if none give does all.",nargs='*')
    cmd_args=argp.parse_args()
    funcs_to_use=cmd_args.funcs
    if "pass" in funcs_to_use:
        quit()



    #required args for functions, lists are handled as OR statements, tuples within lists as AND
    #list of tuples, first element is the list of required arguments and second is the defaults if argument is not found, leaving it as None skips defaults

    required_args ={
        "plot_vdf":[(["coordre","coordinates","cellids"],["coordre=REPLACECOORDRE"]),([("filedir","step")],[None])],
        "plot_vdf_profiles":[(["coordre","coordinates","cellids"],["coordre=REPLACECOORDRE"]),([("filedir","step")],[""])],
        "plot_isosurface":[([("surf_step","surf_var")],["surf_step=10","surf_var='vg_rho'"]),([("filedir","step")],[""])]
        
    }



    #This can be v4 or v5, these both get added into nonrestartcalls and v5nonrestartcalls, here for cleanliness sake
    agnostic_call = [


    # Input and output methods, nooverwrite
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, outputdir=outputLocation+'/'+REPLACEINDEX+'_')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, outputfile=outputLocation+REPLACEINDEX+'_outputfiletest.png', nooverwrite=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, outputfile=outputLocation+REPLACEPREVINDEX+'_outputfiletest.png', nooverwrite=1)",
    "pt.plot.REPLACEFUNC(filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX)",

    # cellids, coordinates
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACECELLID)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACECOORDINATES)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellids=REPLACEMULTIPLECELLID)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordinates=REPLACEMULTIPLECOORDINATES)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, coordre=REPLACEMULTIPLECOORDRE)",



    # Thickness, scale
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=2.0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=2.0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, highres=True)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, scale=2.0, highres=True)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, thick=2.0, highres=True)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, highres=3)",

    # Tick interval
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=0.5,axisunit=6)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, tickinterval=1,axisunit=6)",

    # msec musec titles
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='msec')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='musec')",

    # Watermarks
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmarkb=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NE')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='NW')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SE')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wmark='SW')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, Earth=True)",




    # title, axes, noborders
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title=r'$\mathcal{Title}$ and so forth $\odot$', cbtitle=r'$\mathcal{Color}$')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noborder=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noxlabels=1,noborder=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noborder=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',cbtitle='',noylabels=1,noxlabels=1,noborder=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,nocb=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,internalcb=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, title='',noylabels=1,noxlabels=1,noborder=1,internalcb=1,highres=True)",



    # Overplots and flux lines
    "pt.plot.REPLACEFUNC(filedir=fileLocation, step=REPLACETIME, run=verifydir+REPLACEINDEX)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved=1, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxfile=fluxLocation+fluxname, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, fluxdir=fluxLocation, step=REPLACETIME, fluxthick=0.5, fluxlines=10)",
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, fluxdir=fluxLocation, fluxthick=5, fluxlines=2)",

    # Vscale
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vmin=7.e3, vmax=7.e6)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vmin=7.e-3, vmax=7.e0, vscale=1e-6)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vmin=7.e6, vmax=7.e9, vscale=1e3)",


    # Zoom and units
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50],axisunit=3)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxre=[-10,10,5,50],axisunit=6)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxm=[-10e6,50e6,-5e6,15e6])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, boxm=[-10e6,50e6,-5e6,15e6],axisunit=0)",


    # slicethick (Mostly for vdf,vdf_prof,vdfdiff)
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=0, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=2, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=4, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, slicethick=1e3, coordre=REPLACECOORDRE)",
    # cellsize (same)
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=0.5, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=2, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cellsize=4, coordre=REPLACECOORDRE)",
    # fmin, fmax, setThreshold (same)
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fmax=1.e-12, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fmin=1.e-14,fmax=1.e-12, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-20, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=1.e-15, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, setThreshold=0, coordre=REPLACECOORDRE)",
    # Biglabels(vdf,vdfdiff)
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='A', coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='B', biglabloc=0, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='C', biglabloc=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='D', biglabloc=2, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, biglabel='E', biglabloc=3, coordre=REPLACECOORDRE)",
    # cbulk, center, bvector(same)
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, cbulk=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[-7e5,0,0], coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bvector=1, center=[2e5,2e5,2e5], coordre=REPLACECOORDRE)",

    # wflux(same)
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, wflux=1, coordre=REPLACECOORDRE)",
    # directions(about same)
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xy=1,coordswap=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,0,5], coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, xz=1,coordswap=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[0,1,0], coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, yz=1,coordswap=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[-1,0,0], coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bpara=1,coordswap=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, bperp=1,coordswap=1, coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, normal=[1,1,1], coordre=REPLACECOORDRE)",



    # colormaps 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='nipy_spectral')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='jet')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='hot_desaturated_r')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='viridis')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='plasma')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='magma')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='warhol')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='bwr')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, colormap='PuOr')",

    #colorbar
    #not yet in in master see PR #359
    #"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, cb_horizontal=True)",

    #AMR, fsaved
    #Does not work currently, fix for the AMR contour is in PR #364
    #"pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, amr=0.1,amrlinestyles='dashed',amrcolours='red',amrlinewidths=1"),
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fsaved='red')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxrope=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, fluxrope=0.5)",
    ]


    #This is here so the agnostic calls are first and the order doesnt change if you add something to nonrestartcalls
    nonrestartcalls=[]
    v5nonrestartcalls=[]
    restartcalls=[]
    v5restartcalls=[]
    nonrestartcalls.extend(agnostic_call)
    v5nonrestartcalls.extend(agnostic_call)
    restartcalls.extend(agnostic_call)
    v5restartcalls.extend(agnostic_call)

    v5restartcalls.extend([

    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v')"

    ])

    restartcalls.extend([
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='V')"

    ])

    nonrestartcalls.extend([

    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', colormap='hot_desaturated')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', colormap='hot_desaturated', vscale=1e3)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', op='x', colormap='RdBu',symlog=0, usesci=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', op='y', colormap='RdBu',symlog=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='E', op='z', colormap='RdBu',symlog=0, usesci=0)",

    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='RhoBackstream', colormap='jet')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='beta',lin=1, usesci=0, colormap='viridis',vmax=50)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='MA',lin=1,usesci=0,vmin=2,vmax=40, colormap='inferno_r')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Mms',lin=1,usesci=0, vmin=0, colormap='magma')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='va',lin=1,usesci=0, colormap='magma_r')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vms',lin=1, colormap='plasma_r')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vs',lin=1, colormap='viridis_r')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='max_v_dt', vscale=1e6)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='max_v_dt', vscale=1e3)",

    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Temperature', colormap='plasma', vscale=1e-6,lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pressure', vscale=1e9)",

    # Symlog and vscale
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1e-9)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1e-12)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=0,vscale=1e9)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1,vscale=1e9)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='B', op='y', colormap='bwr',symlog=1e-3,vscale=1e9)",
    
    '''
    # Externals and expressions
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, pass_vars=['rho','B','beta'])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, pass_vars=['va'], vmin=1, vmax=20,lin=1,usesci=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, boxre=[0,30,-15,15], vmin=1, vmax=20,lin=1,usesci=0)",


    # Everything at once
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15], expression=exprMA_cust, vmin=1, vmax=20,lin=1,usesci=0, fsaved=1, fluxfile=fluxLocation+fluxname)",
    '''

    # Streamlines, vectors
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B', vectorcolormap='viridis')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B', vectorcolormap='magma')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=20)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=400)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=20, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='B',vectordensity=400, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectorsize=1.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20,vectorsize=1.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400,vectorsize=1.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',boxre=[-10,10,5,50],vectorsize=1.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20, boxre=[-10,10,5,50],vectorsize=1.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400, boxre=[-10,10,5,50],vectorsize=1.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectorsize=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20,vectorsize=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400,vectorsize=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',boxre=[-10,10,5,50],vectorsize=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=20, boxre=[-10,10,5,50],vectorsize=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='E',vectordensity=400, boxre=[-10,10,5,50],vectorsize=0.5)",


    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinecolor='black')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinecolor='gray')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=0.5, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='B', streamlinedensity=2, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5,streamlinethick=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinethick=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2,streamlinethick=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5,streamlinethick=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinethick=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2,streamlinethick=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=0.5, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='V', streamlinedensity=2, boxre=[-10,10,5,50])",

    # More data reducers
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VBackstream',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VNonBackstream',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VParallel', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VPerpendicular',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VParallelBackstream', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VPerpendicularBackstream',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VParallelNonBackstream', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='VPerpendicularNonBackstream',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pressure')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PParallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpOverPar', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PParallelBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpendicularBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpOverParBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PParallelNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpendicularNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='PPerpOverParNonBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pdyn',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Pdynx',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Temperature')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TParallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TParallelBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpendicularBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TParallelNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpendicularNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpOverPar', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpOverParBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='TPerpOverParNonBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpOverPar', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpOverParBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpOverParNonBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='beta')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaParallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='betaPerpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vBeam',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vBeamRatio')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Thermalvelocity',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='Blocks')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='gyrotropy')"])







    v5nonrestartcalls.extend([
    # Variables, operators, colormaps, usesci, lin, vscale
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', colormap='nipy_spectral')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', colormap='nipy_spectral', vscale=1e9)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', op='x', colormap='bwr')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_b_vol', op='z', colormap='bwr')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', colormap='warhol',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', colormap='warhol',lin=1, vscale=1e-3)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', op='x', colormap='PuOr')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', op='y', colormap='PuOr',symlog=0, usesci=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v', op='z', colormap='PuOr',symlog=0, usesci=0)",

    '''
    # Externals and expressions
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, pass_vars=['vg_rho','vg_b_vol','vg_beta'])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, pass_vars=['vg_va'], vmin=1, vmax=20,lin=1,usesci=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, expression=exprMA_cust, boxre=[0,30,-15,15], vmin=1, vmax=20,lin=1,usesci=0)", #why error with plot3d? keyeerror vg_va
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3, pass_vars=['vg_rho','vg_b_vol','vg_beta'],lin=1,colormap='bwr',usesci=0)",
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=expr_cav_cust, pass_times=3,lin=1,colormap='bwr',usesci=0, boxre=[0,30,-15,15])",
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[7,0], pass_vars=['vg_rho'], boxre=[0,30,-15,15])",
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[7,0], pass_vars=['vg_beta'])",


    # Everything at once
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, external=extcontour, boxre=[0,30,-15,15], expression=exprMA_cust, vmin=1, vmax=20,lin=1,usesci=0, fsaved=1, fluxfile=fluxLocation+fluxname)",
    '''
    
    # Streamlines, vectors
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,vectordensity=200)", 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,normal='x',vectordensity=200)", 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,normal='y',vectordensity=200)", 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectorsize=1,normal='z',vectordensity=200)", 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectordensity=400, boxre=[-10,10,5,50],vectorsize=1.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_v',vectordensity=20, boxre=[-10,10,5,50],vectorsize=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol', vectorcolormap='viridis')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectordensity=400)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectordensity=20, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,vectordensity=200)", 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,normal='x',vectordensity=200)", 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,normal='y',vectordensity=200)", 
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, vectors='vg_b_vol',vectorsize=1,normal='z',vectordensity=200)", 

    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinecolor='black')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinecolor='gray')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=0.5, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_b_vol', streamlinedensity=2, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5,streamlinethick=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinethick=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2,streamlinethick=2)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5,streamlinethick=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinethick=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2,streamlinethick=0.5)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=0.5, boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', boxre=[-10,10,5,50])",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, streamlines='vg_v', streamlinedensity=2, boxre=[-10,10,5,50])",


    # More data reducers
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_rho')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_nonthermal',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_thermal',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_parallel', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_perpendicular',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_parallel_nonthermal', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_perpendicular_nonthermal',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_parallel_thermal', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_v_perpendicular_thermal',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_pressure')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_parallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_perpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_anisotropy', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_parallel_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_perpendicular_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_anisotropy_nonthermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_parallel_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_perpendicular_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_p_anisotropy_thermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_pdyn',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_pdynx',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_temperature')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_parallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_perpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_parallel_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_perpendicular_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_parallel_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_perpendicular_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_anisotropy', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_anisotropy_nonthermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_t_anisotropy_thermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_anisotropy', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_anisotropy_nonthermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_anisotropy_thermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_parallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_beta_perpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='vg_rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",
    "pt.plot.REPLACEFUNC(filename=fileLocation+bulkname, run=verifydir+REPLACEINDEX, expression=timesmooth, pass_times=[14,0],pass_vars=['vg_rho'])"

    ])


    v5multipopcalls=[
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_v_parallel', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_v_perpendicular',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_pressure')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_parallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_perpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_anisotropy', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_parallel_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_perpendicular_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_anisotropy_nonthermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_parallel_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_perpendicular_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_p_anisotropy_thermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_pdyn',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_pdynx',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_temperature')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_parallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_perpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_parallel_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_perpendicular_nonthermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_parallel_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_perpendicular_thermal')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_anisotropy', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_anisotropy_nonthermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_t_anisotropy_thermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_anisotropy', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_anisotropy_nonthermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_anisotropy_thermal', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_parallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_beta_perpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_thermalvelocity',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_blocks')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/vg_gyrotropy')"
    ]


    multipopcalls=[
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, pop='REPLACEPOP', coordre=REPLACECOORDRE)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/VParallel', op='magnitude',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/VPerpendicular',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Pressure')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PParallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpOverPar', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PParallelBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpendicularBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpOverParBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PParallelNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpendicularNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/PPerpOverParNonBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Pdyn',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Pdynx',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Temperature')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TParallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TParallelBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpendicularBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TParallelNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpendicularNonBackstream')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpOverPar', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpOverParBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/TPerpOverParNonBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpOverPar', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpOverParBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpOverParNonBackstream', vmin=0.1, vmax=10)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/beta')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaParallel')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/betaPerpendicular')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Rmirror',lin=1,vmin=0.5,vmax=1.5,usesci=0)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Thermalvelocity',lin=1)",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/Blocks')",
    "pt.plot.REPLACEFUNC(vlsvobj=f, run=verifydir+REPLACEINDEX, var='REPLACEPOP/gyrotropy')"]


    calls = []



    functions = funcs_to_use
    if not functions:
        quit()

    skipped_args=[]



    for j,func in enumerate(functions):
        copy_of_lists=[nonrestartcalls.copy(),restartcalls.copy(),multipopcalls.copy(),v5restartcalls.copy(),v5nonrestartcalls.copy(),v5multipopcalls.copy()]
        for calls_in in copy_of_lists:
            temp_list=[]
            v5restart_replace=False
            restart_replace=False
            if calls_in == v5restartcalls:
                v5restart_replace=True
            elif calls_in == restartcalls:
                restart_replace=True

            for call in calls_in:
                call=call_replace(call,func,skipped_args)
                if not call:
                    continue

                if v5restart_replace:
                    call = call.replace("var='vg_v'","var='vg_restart_v'")
                elif restart_replace:
                    call = call.replace("var='V'","var='restart_V'")

                if call not in calls:
                    calls.append(call)
                    temp_list.append(call)

            calls_in[:]=temp_list



        with open(f"testpackage_{func}.py","w") as f:
                f.write("nonrestartcalls=["+',\n'.join(['"{}"'.format(call) for call in copy_of_lists[0]])+"]\n\n")
                f.write("restartcalls=["+',\n'.join(['"{}"'.format(call) for call in copy_of_lists[1]])+"]\n\n")
                f.write("multipopcalls=["+',\n'.join(['"{}"'.format(call) for call in copy_of_lists[2]])+"]\n\n")

                #v5
                f.write("v5restartcalls=["+',\n'.join(['"{}"'.format(call) for call in copy_of_lists[3]])+"]\n\n")
                f.write("v5nonrestartcalls=["+',\n'.join(['"{}"'.format(call) for call in copy_of_lists[4]])+"]\n\n")
                f.write("v5multipopcalls=["+',\n'.join(['"{}"'.format(call) for call in copy_of_lists[5]])+"]\n\n")


