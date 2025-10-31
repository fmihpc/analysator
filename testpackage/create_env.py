import os
import venv
import subprocess
from sys import version_info as python_version_info




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

def create_venv(path,install_analysator=True,editable=False):
    virt_env= venv.EnvBuilder(with_pip=True,upgrade_deps=True)
    context=virt_env.ensure_directories(path)
    virt_env.create(path)
    virt_env.setup_python(context)

    #Does not work in python versions <3.13
    if python_version_info.major>=3 and python_version_info.minor>=13:
        virt_env.create_git_ignore_file(context)

    virt_env.create_configuration(context)
    virt_env.setup_scripts(context)
    if install_analysator:
        system_call(f'{path}/bin/pip install {'--editable' if editable else ''} ../',live_output=True)
    print(f'Virtual environment created at {path}')
    return None

def gitcheck():

    git_branch = system_call('git rev-parse --abbrev-ref HEAD')


    #Are we on master branch?
    if git_branch not in ('master','main') :
        user_input=input('::warning not in master branch, are you sure you want to continue? y/n\n')
        if user_input.capitalize() != 'Y':
            quit()
    else:
        #Is local master up to date with remote?
        git_diff=system_call('git diff --numstat origin/master...')
        if git_diff!='':
            print(git_diff)
            print('::warning Local master branch not up to date with remote, are you sure you want to continue? y/n\n')
            user_input=input()
            if user_input.capitalize() != 'Y':
                quit()
    return None


'''
def create_ci_venv():
    print('Creating virtual environment for CI')
    import time
    gmtime_now=time.strftime("%c",time.gmtime()).replace(" ","_").replace(":","-")

    gitcheck()

    ci_path="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage"

    path=f"{ci_path}/venv_testpackage-{gmtime_now}-{git_branch}"

    create_venv(path,install_analysator=True,editable=False) #Editable should be false (default) for CI!

    #Make group writable
    print(system_call(f'chmod -R g+w {path}'))

    return path
'''


def create_venv_script(path,venv_path):
    if path[-1]!="/":
        path=path+"/"

    if not os.path.isdir(path):
        raise FileNotFoundError(f"{path} does not exist.")
    
    if os.path.isfile(path+"pyvenv.sh"):
        #If file exists, check for source line and add one if not present
        source_line=False
        with open(path+"pyvenv.sh","a+") as f:
            f.seek(0)
            for line in f:
                if "source" in line[:6]:
                    source_line=True
                    print("Warning: source line already in pyvenv.sh, make sure it is correct for CI! Run regardless? y/n")
                    user_input=input()
                    if user_input.capitalize() != 'Y':
                        quit()
                    break
            if not source_line:
                f.write(f"source {path}/bin/activate\n")
            f.close()
    else:
        #Create the file with the source line
        with open(path+"pyvenv.sh","w") as f:
            f.write("module purge\n")
            if 'TURSO' in os.uname().nodename.upper():
                f.write("export PATH=/wrk-vakka/group/spacephysics/proj/appl/tex-basic/texlive/2023/bin/x86_64-linux:$PATH\n")
            f.write("module load ImageMagick/7.1.0-37-GCCcore-11.3.0\n")
            f.write("module list\n")
            f.write(f"source {path}/bin/activate\n")
            f.close()
    

if __name__ == "__main__":
    
    create_venv_local = True        #Create a virtual environment in the current folder (should be in testpackage folder)


    venv_name= 'venv_testpackage'

    if not venv_name in os.listdir('.') and create_venv_local:
        print('venv_testpackage not found, creating virtual environment')
        create_venv(venv_name,editable=True)
        create_venv_script('./',f'./{venv_name}')
    else:
        print('venv_testpackage found, not creating virtual environment')


    
