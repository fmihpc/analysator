import os
import venv
from sys import version_info as python_version_info


#CHANGE TO NAOTHER FILE
create_venv_CI = False
create_venv_local = True


if not 'venv_testpackage' in os.listdir('.') and create_venv_local:
    print('venv_testpackage not found, creating virtual environment')
    os.system("mkdir venv_testpackage")
    virt_env= venv.EnvBuilder(with_pip=True,upgrade_deps=True)
    context=virt_env.ensure_directories('venv_testpackage')
    virt_env.create('./venv_testpackage')
    virt_env.setup_python(context)

    #Does not work in python versions <3.13
    if python_version_info.major>=3 and python_version_info.minor>=13:
        virt_env.create_git_ignore_file(context)

    virt_env.create_configuration(context)
    virt_env.setup_scripts(context)
    os.system('./venv_testpackage/bin/pip install --editable ../')



