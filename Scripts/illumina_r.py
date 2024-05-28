import sys
import os
import subprocess
import platform


project = sys.argv[1]; workdir = sys.argv[2]; fastq_files = sys.argv[3]; region = sys.argv[4]

if platform.system() == "Windows":

    if len( os.listdir("C:/Program Files/R")) == 0:
        print("The program R is not installed in this computer, please install it first")
        exit()

    else:
        r_script_path = os.path.join("C:/", "Program Files", "R", os.listdir("C:/Program Files/R")[0], "bin", "Rscript")
        subprocess.call([r_script_path, "Scripts/DADA_illumina_V3V4_V4.R", project, workdir, fastq_files, region], shell=False)

elif platform.system() == "Linux":
    subprocess.call(["Rscript", "Scripts/DADA_illumina_V3V4_V4.R",project, workdir, fastq_files, region], shell=True)

else:
    print(f"Platform {platform.system()} not supported")
    exit()


