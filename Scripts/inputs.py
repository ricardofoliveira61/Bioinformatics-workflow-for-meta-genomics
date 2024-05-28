import re; import os; import platform; import subprocess; import os


def action():
    inp = int(input("What do you pretend to do [1:amplicon analysis, 2: ASV classification, 3: phyloseq analysis, 4: Exit](pick a number): "))
    if inp in range(1,5): return inp
    else: inp = int(input("Invalid answer. Choose from 1:amplicon analysis, 2: ASV classification, 3: phyloseq analysis, 4: Exit "))


def script(state:int,workdir):
    
    if state == 1: 
        project, fastq_files, software, seq_tech, region = get_user_inputs ()
        amplicon_script_2_run(project, workdir, fastq_files, software, seq_tech, region)
    
    elif state == 2:
        project, classifier_met, data_base, fasta_path = classifier()
        classier_scrit_2_run(workdir, project, classifier_met, data_base, fasta_path)
        

    elif state == 3:
        print("Not available yet")

    else:
        exit()


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
    project = input("How do you want to name your files?: ")
    while True:
        if not re.match(r'^[0-9A-Za-z_-]+$',project):
            project = input("Invalid file name. Use only alphanumeric characters (a-z, A-Z, 0-9), underscores (_) and hyphens (-): ")
        
        else: return project


def work_dir():
    workdir = input("What is the working directory? (Must be the absolute path): ")
    while True:
        if os.path.isabs(workdir): return workdir
        else: workdir = input ("Invalid working directory: ")



def direcoty_paths ():
    fastq_files = input("What is the diretory of the fastqfiles? (Must be the absolute path): ")
    while True:
        if os.path.isabs(fastq_files): return fastq_files
        else: fastq_files = input ("Invalid directory: ")        


def analysis_software():

    software = input("Which software you want to use? [DADA2, Vsearch]: ").upper()

    while True:
        if software in ["DADA2", "VSEARCH"]: return software
        else: software = input ("Invalid software. Please choose from: DADA2 or Vsearch: ").upper()


def get_user_inputs ():

    project = project_name()
    fastq_files = direcoty_paths() 
    software = analysis_software()
    seq_tech, region = region_seqtech_input()
    return project, fastq_files, software, seq_tech, region


def dada2_illumina(project:str, workdir:str, fastq_files:str, region:str):
    

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


def amplicon_script_2_run ( project:str, workdir:str, fastq_files:str, sofware:str, seq_tech:str, region:str):
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


# def what_next():
#     ans = int(input("\n\nWhat do want to do next [1: New amplicon analysis, 2: Procced to amplicon classification, 3: Exit](Use only numbers): "))
#     while True:
#         if ans in [1, 2, 3]: return ans
#         else: ans = int(input ("Invalid answer. Please choose the number from: 1: New amplicon analysis, 2: Procced to amplicon classification, 3: Exit: "))


def classifier_method():
    classifier = input("Which classifier you intend to use [Decipher, rdp, sintax]: ").upper()
    while True:
        if classifier in ["DECIPHE", "RDP","SINTAX"]: return classifier
        else: classifier = input("Invalid classifier. Please choose from Decipher, rdp or sintax: ").upper()


def database():
    data_base = input("Which database you want to use [Silva, RDP, greengenes]: ").upper()
    while True:
        if data_base in ["SILVA", "RDP", "GREENGENES"]: return data_base
        else: data_base = input("Invalid database. Please choose from Silva, RDP or greengenes: ").upper()


def fasta_file():
    fasta_path = input("What is the absolut path to the ASV fasta file (must include the fasta file): ")
    while True:
        if os.path.isabs(fasta_path) and fasta_path.endswith((".fasta", ".fa",".fas",".fsa")): return fasta_path
        else: fasta_path = input ("Invalid directory: ")


def classifier():
    project = project_name()
    classifier_met = classifier_method()
    data_base = database()
    fasta_path = fasta_file()
    return project, classifier_met, data_base, fasta_path


def classier_scrit_2_run(workdir, project, classifier_met, data_base, fasta_path):
    if classifier_met == "DECIPHER":
        if platform.system() == "Windows":

            if len( os.listdir("C:/Program Files/R")) == 0:
                print("The program R is not installed in this computer, please install it first")
                exit()

            else:
                r_script_path = os.path.join("C:/", "Program Files", "R", os.listdir("C:/Program Files/R")[0], "bin", "Rscript")
                print("\n\nStarting ASVs assignment using DECIPHER")
                subprocess.call([r_script_path, "Scripts/decipher_rdp_classifiers.R", project, workdir, data_base, fasta_path], shell=False)
                print("\n\nTask concluded")

        elif platform.system() == "Linux":
            print("\n\nStarting ASVs assignment using DECIPHER")
            subprocess.call(["Rscript", "Scripts/decipher_rdp_classifiers.R",project, workdir, data_base, fasta_path], shell=True)
            print("\n\nTask concluded")

    else:
        print(f"Platform {platform.system()} not supported")
        exit()



#if __name__ =="__main__":
#    script(1,"C:\\Users\\ricar\\Desktop\\Teste")

