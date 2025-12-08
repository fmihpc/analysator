import subprocess
import time as timetime
import os
import sys
from argparse import ArgumentParser
import cv2
import numpy as np
import logging
import os.path

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



def compare_images_in_folders(a,b,jobcount,jobcurr):


    #folders=get_directories(a,depth=2)
    #folders2=get_directories(b,depth=2)
    
#    if len(folders)!=len(folders2):
#        raise SystemError("Folders do not have the same structure!")

    #do the comparisons
    
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


    nteststot = len(different_files)
    increment = int(nteststot/jobcount)
    remainder = nteststot - jobcount * increment
    start=jobcurr * increment
    end=start + increment
    # Remainder frames are divvied out evenly among tasks
    if jobcurr < remainder:
        start = start + jobcurr
        end = end + jobcurr + 1
    else:
        start = start + remainder
        end = end + remainder



    #Feed the different files to compare_images
    for file in different_files[start:end]:
        if(not compare_images(file,file.replace(a,b))):
            different=True
            print("Images differ:",file,file.replace(a,b))



    #Print unique and missing files
    for file in unique_files:
        print("Unique file:",file)
    for file in missing_files:
        print("Missing file:",file)


    if len(unique_files)!=0:
        print("::warning::Found new file(s) produced by the code!")

    if len(missing_files)!=0:
        raise SystemError("Found file(s) **not** produced by the code!")


    if different:
        print(f"::error title=Plot(s) differ::Produced plots not in agreement with the verfication set {a}")
        raise SystemError("Images Differ")


#only prints what is at the target depth (might be clunky as i was tired when writing it)
def get_directories(directory,depth):
    out_dirs=[]
    dirs=[]
    if depth==0:
        return [directory]
    for i in range(depth):

        if i==0:
            directory= [directory]
            subdirs_current_depth=os.listdir(directory[0])
            for subdir in subdirs_current_depth:
                dirs.append(os.path.join(directory[0],subdir))
            directory=dirs
        else:
            directory=[]

            for dircs in dirs:
                subdirs_current_depth=os.listdir(dircs)
                for subdir in subdirs_current_depth:
                    directory.append(os.path.join(dircs,subdir))
            length_current_depth=len(directory)
            dirs=directory
   
        if i==depth-1:
            out_dirs=directory

    return out_dirs


if __name__=='__main__':

    #Parse arguments
    parser = ArgumentParser(prog="Image compare"
                            ,description="Compares images in two folders/images, and outputs the different ones to a specified output folder. Note that folders must have the same structure and filenames, otherwise these files are treated as unique."
                            )

    parser.add_argument("folder_a",help="First folder/image to compare")
    parser.add_argument("folder_b",help="Second folder/image to compare")
   
    parser.add_argument("jobcount",help="Number of parallel jobs to use",default=1,nargs='?',type=int)
    parser.add_argument("jobindex",help="Index of the job to run",default=0,nargs='?',type=int)
    
   
    args= parser.parse_args()
    a,b = args.folder_a,args.folder_b

    jobcount=args.jobcount
    jobindex=args.jobindex

    compare_images_in_folders(a,b,jobcount,jobindex)