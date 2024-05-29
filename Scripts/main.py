import inputs; import os



workdir = inputs.work_dir()
main_path = os.path.dirname(os.path.abspath(__file__))
while True:

    state = inputs.action()
    inputs.script(state, workdir, main_path)
    print()






    # if first_run:
    #     project, workdir, fastq_files, sofware, seq_tech, region = inputs.get_user_inputs(first_run); first_run = False

    # else: project, fastq_files, sofware, seq_tech, region = inputs.get_user_inputs(first_run)   

    # inputs.script_2_run(project, workdir, fastq_files, sofware, seq_tech, region) 

    # ans =inputs.what_next()
    # if ans == 1:
    #     continue
    # elif ans == 3:
    #     exit()

    # print("Classification not avaivle for now")
    # exit()
