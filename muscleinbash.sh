#!/bin/bash 

echo Running MUSCLE file...

chmod +x MUSCLE_HL.R
Rscript --slave --vanilla --no-restore MUSCLE_HL.R

