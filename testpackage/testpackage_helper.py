import subprocess

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

