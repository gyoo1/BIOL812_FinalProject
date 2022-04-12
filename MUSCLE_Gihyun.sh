{\rtf1\ansi\ansicpg1252\cocoartf2638
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fnil\fcharset0 HelveticaNeue;\f1\fnil\fcharset0 HelveticaNeue-Bold;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab560
\pard\pardeftab560\slleading20\partightenfactor0

\f0\fs26 \cf0 #!/bin/bash\
#SBATCH -c 4                               # Number of CPUS requested. If omitted, the default is 1 CPU.\
#SBATCH --mem=2000                    # Memory requested in megabytes. If omitted, the default is 1024 MB.\
#SBATCH -t 2:0:0      # How long will your job run for? If omitted, the default is 3 hours.\
	\
module load StdEnv/2020 r/4.1.2	\
\
\pard\pardeftab560\slleading20\partightenfactor0

\f1\b \cf0 echo
\f0\b0  Running MUSCLE file...\
\

\f1\b chmod
\f0\b0  +x MUSCLE_Gihyun.R\
Rscript
\f1\b  --slave --vanilla --no-restore
\f0\b0  MUSCLE_Gihyun.R\
}