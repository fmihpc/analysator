import os
import venv
import subprocess
from sys import version_info as python_version_info


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
        os.system(f'{path}/bin/pip install {'--editable' if editable else ''} ../')
    print(f'Virtual environment created at {path}')
    return 0

def system_call(cmd,wait=False):
    proc = subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    out,err = proc.communicate()

    #If errors, raise an exception
    if err:
        err = str(err,'utf-8')
        raise RuntimeError(err)

    out = str(out,'utf-8').rstrip('\n')
    if wait:
        proc.wait()
    return out


if __name__ == "__main__":
    create_venv_CI = True
    create_venv_local = False

    ci_path="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage"

    if not 'venv_testpackage' in os.listdir('.') and create_venv_local:
        print('venv_testpackage not found, creating virtual environment')
        os.system("mkdir venv_testpackage")  #Is this necessary?
        create_venv('venv_testpackage',editable=True)

    if create_venv_CI:
        print('Creating virtual environment for CI')
        import time
        time=time.strftime("%c",time.gmtime()).replace(" ","_").replace(":","-")

        git_branch = system_call('git rev-parse --abbrev-ref HEAD')


        #Are we on master branch?
        if git_branch!='master':
            user_input=input('::warning not in master branch, are you sure you want to continue? y/n\n')
            if user_input.capitalize() != 'Y':
                quit()

        path=f"{ci_path}/venv_testpackage-{time}-{git_branch}"

        create_venv(path,install_analysator=True,editable=False) #Editable should be false (default) for CI!
        
        #Make group writable
        print(system_call(f'chmod g+w {path}'))
