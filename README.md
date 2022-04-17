
#BIOL812 Final Group Project
##Authors: Gihyun Yoo, Heather Lounder, Kayleigh Casmey, Glafira Ermakova

This is an R Markdown document that generates a github readme.md file.

Assignment prepared for BIOL 812 Introduction to Computational Analysis Winter 2022. In this pipeline, we concatenated the DNA sequences for several Coronavirus variants (Including MERS, SARS and COVID-19) in order to run a MUSCLE multiple sequence alignment. Using this alignment we generated a distance matrix and a phylogeny tree in order to examine the evolutionary relationship between the different variants. 

#Scripts
1. combining_seq.R/py - converts the chosen FASTA sequences into 1 concatenated DNA string file needed for MUSCLE
2. MUSCLE.R - MUSCLE sequence alignment script for the concatenated DNA string file
3. MUSCLEUNIX.sh - runs the MUSCLE.R script in commandline 
4. Tree.R - Creates a distance matrix and annotated phylogeny tree based on the MUSCLE sequence alignment 

![pipeline for the BIOL 812 assignment](./Pipeline.png)

