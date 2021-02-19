#!/bin/tcsh
#This file is called submit-script.sh
#  #SBATCH --partition=centos7
#SBATCH --time=4:00:00   # run time in days-hh:mm:ss
#SBATCH --ntasks=16          # default 16 if this line not specified
#SBATCH --mem=128gb
##SBATCH --mem-per-cpu=2gb
#SBATCH -o out.slurm-%j.%N #
#SBATCH -e err.slurm-%j.%N #
#Make sure to change the above two lines to reflect your appropriate
# file locations for standard error and output

source /p/stellopt/ANALYSIS/jschmitt/src/module_fotd_regcoil_vc7

mpirun -n 16 --mca mtl ^psm --mca btl ^openib  /p/stellopt/ANALYSIS/jschmitt/src/STELLOPT_REGCOIL/STELLOPTV2/Release/xstelloptv2 input.REGCOIL_UNIFORM >& log.stellopt



