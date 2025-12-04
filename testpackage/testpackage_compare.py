import subprocess
import time as timetime
import os
import sys
from argparse import ArgumentParser
import cv2
import numpy as np
import logging

#could be used to replace compare_images with a cv2 based implementation
def compare_images(a,b):

    im1=cv2.imread(a).astype(np.float16)
    im2=cv2.imread(b).astype(np.float16)
    #originally uint8 so we might underflow with the substraction if not casted as float, 16 should be good enough
    if im1.shape != im2.shape:
        '''
        something should be added to handle this better, have a threshold in general or something.
        the substraction yields an error, so one could pad it with 
        smaller_image=np.pad(smaller_image,(0,larger_image.shape[0]-smaller_image.shape[0]),(0,larger_image.shape[1]-smaller_image.shape[1]),(0,0))
        '''
        return False
    diff = np.sqrt(np.mean((im1-im2)**2))
    if diff !=0:
        return False
    return True



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
    missing_files=[]
    out=out.splitlines()

    #Parse the output of diff
    for line in out:
        line=line.split(" ")
        if line[0]=="Binary" and line[1]=="files":
            different_files.append(line[2])
        elif line[0]=="Only":
            if b in line[2]:
                unique_files.append(line[2].rstrip(":")+'/'+line[3])
            else:
                missing_files.append(line[2].rstrip(":")+'/'+line[3])


    #Feed the different files to compare_images
    for file in different_files:
        #if output_folder!= "NULL:":
        #    filename = file.split("/")[-1].rstrip(".png") #is it always png?
        #    output_folder=output_folder+f"/difference_output_{filename}.png"

        if(not compare_images(file,file.replace(a,b))):
            different=True
            print("Images differ:",file,file.replace(a,b))



    #Print unique and missing files
    for file in unique_files:
        print("Unique file:",file)
    for file in missing_files:
        print("Missing file:",file)


    if len(unique_files)!=0:
        print("::warning Found new file(s) produced by the code!")

    if len(missing_files)!=0:
        raise SystemError("Found file(s) **not** produced by the code!")


    if different:
        print(f"::error title=Plot(s) differ::Produced plots not in agreement with the verfication set {a}")
        raise SystemError("Images Differ")

if __name__=='__main__':

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

    compare_images_in_folders(a,b,output_folder)
