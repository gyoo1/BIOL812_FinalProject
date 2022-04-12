<<<<<<< HEAD
#Combining FASTA files into single file 

=======
>>>>>>> 3feb83da8228dfb076aeeebc39315bd819f32533
import os
DIR = 'Sequences'
os.chdir(DIR)
DIR = os.getcwd()
file_name = input( "What you you like to name output file? Note: please include .fasta extension: ")

output_fasta= open(file_name, "w")
file_count=0

for f in os.listdir(DIR):
    if f.endswith((".fasta")):
        file_count+=1
        fh= open(os.path.join(DIR, f))
        for line in fh:
            output_fasta.write(line)
        fh.close()

output_fasta.close()
print(str(file_count)) + "FASTA files have been merged into single file, can be found at " + DIR

