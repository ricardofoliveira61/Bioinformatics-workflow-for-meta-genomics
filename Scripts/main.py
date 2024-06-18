import inputs; import os


print("""
      Welcome to the toolbox!
      This toolbox was developed by Ricardo Oliveira and Pedro Santos. 
      For more information or help, please contact us using the follwing email: ricardofoliveira61@gmail.com.

      This toolbox incorporate various famous tools used in metabarcoding pipelines.
      The amplicon analysis is currently performed using the R package DADA2 for more information please visit: https://benjjneb.github.io/dada2/index.
      The taxonomic assignment is done using the R package DECIPHER for more information please visit: http://www2.decipher.codes/
      The results exploration is performed using the R package phyloseq for more information please visit: https://joey711.github.io/phyloseq/index.html

      Problems, errors or suggestions must be reported in the toolbox github through the creation of new issue or pull request repository: https://github.com/ricardofoliveira61/

      Have a good time using the toolbox!
      """)

workdir = inputs.work_dir()
main_path = os.path.dirname(os.path.abspath(__file__))
while True:

    state = inputs.action()
    inputs.script(state, workdir, main_path)
    print()
