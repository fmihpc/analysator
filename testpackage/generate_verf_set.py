import subprocess
import os



def system_call(cmd):
    proc = subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)

    out,err = proc.communicate()

    #If errors, raise an exception
    if err:
        err = str(err,'utf-8')
        raise RuntimeError(err)

    out = str(out,'utf-8').rstrip('\n')
    return out

#output_dir named after the latest commit hash
output_dir="/wrk-vakka/turso/group/spacephysics/CI_analysator/analysator_testpackage/verification_sets/"+system_call("git rev-parse HEAD")


git_branch = system_call('git rev-parse --abbrev-ref HEAD')

#Are we on master branch?
if git_branch!='master':
    user_input=input('::warning not in master branch, are you sure you want to continue? y/n\n')
    if user_input.capitalize() != 'Y':
        quit()

#Check if folder exists
if not os.path.isdir(output_dir):
    system_call(f'mkdir -p {output_dir}')
else:
    print('::warning output folder already exists, we might be overwriting something')

#Call the sbatch
out=system_call(f'sbatch ./run_testpackage_generate_verf_set.sh {output_dir}')
print(out)