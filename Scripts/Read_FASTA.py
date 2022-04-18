#Combining FASTA files into single file 

import os
os.getcwd() #current directory
DIR = input("\nInput folder containing FASTA files to combine : ")
os.chdir(DIR)
FILE_NAME = input("Enter output file name with .fasta extension : ")
output_fas = open(FILE_NAME, 'w')
file_count = 0

for f in os.listdir():
        if f.endswith((".fasta")):
                file_count += 1
                with open(f) as fh:
                        output_fas.writelines(fh)
        fh.close()
print(str(file_count)) + " FASTA files were merged into one file"
