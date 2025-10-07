import subprocess
import time as timetime
import os
import sys
from argparse import ArgumentParser



def compare_images(a,b,output_file="NULL:"):

    cmd = f"magick compare -metric RMSE {a} {b} {output_file}"
    proc = subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.STDOUT)

    out,err = proc.communicate()

    exitcode=proc.returncode

    #If errors, raise an exception (This has to be odne like this as compare sends output to stderr)
    if exitcode!=0 and exitcode != 1:
        out = str(out,'utf-8')
        raise RuntimeError(out)

    out = str(out,'utf-8') 
    out=out.strip('\n')
    out = out.split(" ")

    #Returns true if images' RMSE = 0, i.e are identical
    if out[0]=='0':
        return True

    return False


#output_folder = "/home/siclasse/analysator/different_output/"

#!!!!!outputs unique files only from folder "b"


#Parse arguments
parser = ArgumentParser(prog="Image compare"
                        ,description="Compares images in two folders/images, and outputs the different ones to a specified output folder. Note that folders must have the same structure and filenames, otherwise these files are treated as unique."
                        )

parser.add_argument("folder_a",help="First folder/image to compare")
parser.add_argument("folder_b",help="Second folder/image to compare")
parser.add_argument("output_folder",help="Output folder for different images, if not specified, no output is saved",default="NULL:",nargs='?')

args= parser.parse_args()
a,b,output_folder = args.folder_a,args.folder_b,args.output_folder

#Create output folder if it doesn't exist
if not os.path.exists(output_folder) and output_folder!='NULL:':
    proc = subprocess.Popen(f'mkdir {output_folder}'.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = proc.communicate()

    #If errors, raise an exception
    if err:
        err = str(err,'utf-8')
        raise RuntimeError(err)



def compare_images_in_folders(a,b,output_folder='NULL:'):
    cmd = f'diff -r {a} {b}'
    proc = subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out,err = proc.communicate()
    different=False
    #If errors, raise an exception
    if err:
        err = str(err,'utf-8')
        raise RuntimeError(err)

    out = str(out,'utf-8') 

    different_files=[]
    unique_files=[]
    out=out.splitlines()

    #Parse the output of diff
    for line in out:
        line=line.split(" ")
        if line[0]=="Binary" and line[1]=="files":
            different_files.append(line[2])
        elif line[0]=="Only" and line[2].rstrip(":")==b:
            unique_files.append(line[2].rstrip(":")+line[3])

    #Feed the different files to compare_images
    for file in different_files:
        if(not compare_images(file,file.replace(a,b))):
            different=True

            if output_folder!= "NULL:":
                filename = file.split("/")[-1].rstrip(".png") #is it always png?
                output_folder=output_folder+f"/difference_output_{filename}.png"

            compare_images(file,file.replace(a,b),output_folder)

            print("Images differ:",file,file.replace(a,b))



    #Print unique files
    for file in unique_files:
        print("Unique file:",file)

    if len(unique_files)!=0:
        print("::warning title=Unique file(s)::Found new file(s) produced by the code!")
    if different:
        print(f"::error title=Plot(s) differ::Produced plots not in agreement with the verfication set {a}")
        raise SystemError("Images Differ")


compare_images_in_folders(a,b,output_folder)
