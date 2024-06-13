import inputs; import os


print("Introduction text here. Refer which tools are used")
workdir = inputs.work_dir()
main_path = os.path.dirname(os.path.abspath(__file__))
while True:

    state = inputs.action()
    inputs.script(state, workdir, main_path)
    print()
