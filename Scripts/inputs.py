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
        else: region = input(f"'{region}' is an invalid input. Please choose from: V3 or V3-V4: ").upper()
    
    return seq_tech, region


def project_name():
    project = input("How do you want to name your analysis?: ")
    while True:
        if not re.match(r'^[A-Za-z_-]+$',project):
            project = input("Invalid project name. Use only alphanumeric characters (a-z, A-Z, 0-9), underscores (_) and hyphens (-).: ")
        
        else: return project


def direcoty_paths ():
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
        pass