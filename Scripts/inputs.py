import re
import os


def region_seqtech_input ():

    seq_tech = input("What is the sequence thecnology used? [Illumina, PacBio, ONP]: ").upper()
    while True:
        if seq_tech in ["ILLUMINA", "PACBIO", "ONP"]: break
        else: seq_tech= input(f"'{seq_tech}' is an invalid input. Please choose from: Illumina, PacBio or ONP: ").upper()
    
    region = input("Which region was sequencided? [V4, V3-V4]: " ).upper()
    while True:
        if region in ["V4", "V3-V4"]: break
        else: region = input(f"'{region}' is an invalid input. Please choose from: V4 or V3-V4: ").upper()
    
    return seq_tech, region


def project_name():
    project = input("How do you want to name your analysis?: ")
    while True:
        if not re.match(r'^[A-Za-z_-]+$',project):
            project = input("Invalid project name. Use only alphanumeric characters (a-z, A-Z, 0-9), underscores (_) and hyphens (-): ")
        
        else: return project


def direcoty_paths (first=True):
    if first:
        workdir = input("What is the working directory? (Must be the absolute path): ")
        while True:
            if os.path.isabs(workdir):
                break
            else: workdir = input ("Invalid working directory: ")

        fastq_files = input("What is the diretory of the fastqfiles? (Must be the absolute path): ")
        while True:
            if os.path.isabs(fastq_files):
                break
            else: fastq_files = input ("Invalid directory: ")

        return workdir, fastq_files

    else:
        fastq_files = input("What is the diretory of the fastqfiles? (Must be the absolute path): ")
        while True:
            if os.path.isabs(fastq_files):
                break
            else: fastq_files = input ("Invalid directory: ")

        return fastq_files


def analysis_software():
    software = input("Which software you want to use? [DADA2, Vsearch]: ").upper()

    while True:
        if software in ["DADA2", "VSEARCH"]:
            return software
        else: software = input ("Invalid software. Please choose from: DADA2 or Vsearch: ").upper()


def get_user_inputs (first = True):
    if first:
        project = project_name()
        workdir, fastq_files = direcoty_paths() 
        software = analysis_software()
        seq_tech, region = region_seqtech_input()

        return project, workdir, fastq_files, software, seq_tech, region

    else:
        project = project_name()
        fastq_files = direcoty_paths(first) 
        software = analysis_software()
        seq_tech, region = region_seqtech_input()
        return project, fastq_files, software, seq_tech, region


def dada2_illumina(project:str, workdir:str, fastq_files:str, region:str):
    import platform; import subprocess; import os

    if platform.system() == "Windows":

        if len( os.listdir("C:/Program Files/R")) == 0:
            print("The program R is not installed in this computer, please install it first")
            exit()

        else:
            r_script_path = os.path.join("C:/", "Program Files", "R", os.listdir("C:/Program Files/R")[0], "bin", "Rscript")
            print("\n\nStarting default amplicon analysis using DADA2 for Illumina platform")
            subprocess.call([r_script_path, "Scripts/DADA_illumina_V3V4_V4.R", project, workdir, fastq_files, region], shell=False)

    elif platform.system() == "Linux":
        print("\n\nStarting default amplicon analysis using DADA2 for Illumina platform")
        subprocess.call(["Rscript", "Scripts/DADA_illumina_V3V4_V4.R",project, workdir, fastq_files, region], shell=True)

    else:
        print(f"Platform {platform.system()} not supported")
        exit()


def script_2_run ( project:str, workdir:str, fastq_files:str, sofware:str, seq_tech:str, region:str):
    import sys
    
    if sofware == "DADA2":
        if seq_tech == "ILLUMINA":
            dada2_illumina(project, workdir, fastq_files, region)
            print("\n\nAnalysis concluded. Verify the results before proceeding")
        
        elif seq_tech == "PACBIO":
            pass

        elif seq_tech == "ONT":
            pass

        else:
            print(f"The platform {seq_tech} is not supported")