import subprocess
import analysator as pt #this import is used, see the function_pars
import inspect
import re

def system_call(cmd,live_output=False):
    with subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE) as proc:
        if live_output:
            for line in proc.stdout:
                print(str(line,'utf-8').rstrip('\n')) #Note that for example pip's progress bar is not displayed

        out,err = proc.communicate()
    

    #If errors, raise an exception
    if proc.returncode!=0:
        err = str(err,'utf-8')
        raise RuntimeError(err)

    out = str(out,'utf-8').rstrip('\n')
    return out

def call_replace(call,func,skipped_args,required_args):
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
                if type(param) is tuple and param[0] not in named_parameters:
                    check=True
                elif any((all(r in named_parameters for r in param),(param in named_parameters))):
                    if not skipped_args or func not in skipped_args or not any(r in skipped_args[func].keys() for r in param):
                        check=True
                        break
            if not check:
                #Add parameters if there are default_params
                if default_params:
                    if type(default_params) is str:
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
                    if type(skipped_args_dict[call_args[0]]) is str and skipped_args_dict[call_args[0]] in call_args[1]:
                        continue
                    elif type(skipped_args_dict[call_args[0]]) is list and any(arg_skip in call_args[1] for arg_skip in skipped_args_dict[call_args[0]]):  #list of args in dict value means OR  ex. {'var':["vg_rho","vg_phi"]}
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



