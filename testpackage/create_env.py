import os
import venv
from testpackage_helpers import system_call
from sys import version_info as python_version_info
from sys import version as python_version
import argparse

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
    virt_env.post_setup(context)
    if install_analysator:
        editable='--editable' if editable else None
        system_call(f'{path}/bin/pip install {editable} ../[testpackage]',live_output=True)
    print(f'Virtual environment created at {path}')
    return None


def create_venv_script(path,venv_path):
    if path[-1]!="/":
        path=path+"/"

    if not os.path.isdir(path):
        raise FileNotFoundError(f"{path} does not exist.")
    
    if os.path.isfile(path+"pyvenv.sh"):
        #If file exists, check for source line and add one if not present
        with open(path+"pyvenv.sh","a+") as f:
            f.seek(0)
            for line in f:
                if "source" in line[:6]:
                    raise SystemError("source line already in pyvenv.sh!")

            f.write(f"source {venv_path}/bin/activate\n")
            f.close()
    else:
        #Create the file with the source line
        with open(path+"pyvenv.sh","w") as f:
            f.write("module purge\n")
            if 'TURSO' in os.uname().nodename.upper():
                f.write("export PATH=/wrk-vakka/group/spacephysics/proj/appl/tex-basic/texlive/2023/bin/x86_64-linux:$PATH\n")

            #Get used python version and gcc version
            #note that this may break if the version string format changes
            #Module load of python is required, otherwise python cannot be called outside the directory the venv is in (if the module system is used)
            version_info = python_version.split(" ")
            used_python_version = version_info[0]
            used_gcc_version = version_info[-1].strip("[]")
            
            #If module system is used, on HILE use cray-python!
            if 'HILE' in os.uname().nodename.upper():
                f.write("module load cray-python\n")
            else:
                f.write(f"module load Python/{used_python_version}-GCCcore-{used_gcc_version}\n")
            f.write("module list\n")
            f.write(f"source {venv_path}/bin/activate\n")

            f.close()
    

if __name__ == "__main__":
    #Will install python venv with same version as the python this script was called with 

    create_venv_local = True        #Create a virtual environment in the current folder (should be in testpackage folder)

    parser=argparse.ArgumentParser(description='Create python virtual environment for testpackage, will also create batch script that can be sourced.')
    parser.add_argument('--no-analysator',action='store_true',help='Do not install analysator.',default=False)
    parser.add_argument('--editable','-e',action='store_true',help='Install analysator as editable',default=True)
    args=parser.parse_args()
    
    venv_name= 'venv_testpackage'
    venv_path = os.path.abspath('./'+venv_name)
    if venv_name not in os.listdir('.') and create_venv_local:
        print('venv_testpackage not found, creating virtual environment')
        create_venv(venv_path,editable=args.editable,install_analysator=not args.no_analysator)
        create_venv_script('./',venv_path)
    else:
        print('venv_testpackage found, not creating virtual environment')


    
