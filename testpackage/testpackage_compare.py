import subprocess
import time as timetime
import os
from argparse import ArgumentParser



def compare_images(a,b):
    cmd = f'magick compare -metric RMSE {a} {b} NULL:'
    proc = subprocess.Popen(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    out,err = proc.communicate()

    #If errors, raise an exception
    if err:
        err = str(err,'utf-8')
        raise RuntimeError(err)

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
            filename = file.split("/")[-1].rstrip(".png") #is it always png?
            cmd = f"compare -metric RMSE {file} {file.replace(a,b)} {output_folder}/difference_output_{filename}.png"

            if output_folder=='NULL:':
                cmd = cmd.replace(f"{output_folder}/difference_output_{filename}.png","NULL:")

            proc = subprocess.run(cmd.split(" "),stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            print("Images differ:",file,file.replace(a,b))

            # The error handling doesnt work, the compare sends the result to stderr also, so that's just great
            '''
            out,err = proc.stdout,proc.stderr
            print(out,err)

            #If errors, raise an exception
            if err:
                err = str(err,'utf-8')
                raise RuntimeError(err)
            '''


    #Print unique files
    for file in unique_files:
        print("Unique file:",file)



compare_images_in_folders(a,b,output_folder)
