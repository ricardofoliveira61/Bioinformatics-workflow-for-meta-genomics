import inputs


first_run = True
while True:
    if first_run:
        project, workdir, fastq_files, sofware, seq_tech, region = inputs.get_user_inputs(first_run); first_run = False

    else: project, fastq_files, sofware, seq_tech, region = inputs.get_user_inputs(first_run)   

    inputs.script_2_run(project, workdir, fastq_files, sofware, seq_tech, region) 

    ans =inputs.what_next()
    if ans == 1:
        continue
    elif ans == 3:
        exit()

    print("Classification not avaivle for now")
    exit()
