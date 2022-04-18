#!/bin/bash
#SBATCH -c 4                               # Number of CPUS requested. If omitted, the default is 1 CPU.
#SBATCH --mem=2000                    # Memory requested in megabytes. If omitted, the default is 1024 MB.
#SBATCH -t 2:0:0      # How long will your job run for? If omitted, the default is 3 hours.

module load StdEnv/2020 r/4.1.2

#run MUSCLE

R CMD BATCH MUSCLE.R