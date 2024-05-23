import inputs
import subprocess

first_run = True
while True:
    if first_run:
        project, workdir, fastq_files, sofware, seq_tech, region = inputs.get_user_inputs()        
        first_run = False

    if sofware == "DADA2":

        if seq_tech == "ILLUMINA":
            subprocess.run(["python", "./scripts/illumina_r.py" , project,workdir ,fastq_files, region])


    break








