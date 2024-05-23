import importlib.util
import pip
import sys
import os

if importlib.util.find_spec("rpy2") is not None: print("The package rpy2 is installed.")
else: 
    pip_ = input("The package rpy2 is not installed. Do you want to install it?[Y/N]: ")

    if pip_.upper() == "Y" :
        pip.main(["install","rpy2"])
        print("Package rpy2 installed successfully")
        import rpy2
        from rpy2.robjects import r
    
    else:
        print("The package rpy2 is necessary to carry on amplicon analysis using DADA2")
        exit()


project = sys.argv[1]
workdir = sys.argv[2]
fastq_files = sys.argv[3]
region = sys.argv[4]


# Use the variables in your script
print(f"Received project: {project}")
print(f"Received working directory: {workdir}")
print(f"Received fastq files directory: {fastq_files}")
print(f"Received sequenced region: {region}")


r(f"setwd({workdir})") 