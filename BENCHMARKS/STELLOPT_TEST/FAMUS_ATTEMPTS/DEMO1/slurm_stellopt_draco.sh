#!/bin/bash -l
# Standard output and error:
#SBATCH -o ./tjob.out.%j
#SBATCH -e ./tjob.err.%j
# Initial working directory:
#SBATCH -D ./
# Job Name:
#SBATCH -J j_st
# Queue (Partition):
#SBATCH --partition=express
# Number of nodes and MPI tasks per node:
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=32
#
#SBATCH --mail-type=none
#SBATCH --mail-user=<userid>@rzg.mpg.de
#
#SBATCH --mem=120000
#
# Wall clock limit:
#SBATCH --time=0:30:00

# Load modules and run the program:
#module load hdf5-mpi
#module load netcdf-mpi
#module list
source ~/prep_stellopt_famus.sh
#
# Added for a quick fix. remove when necc.  JCS 2019-01-08    
#export NETCDF_HOME=/mpcdf/soft/SLE_12_SP3/packages/haswell/netcdf-mpi/intel_18_0_3-impi_2018_3/4.4.1/
srun ~/src/STELLOPT_FAMUS/STELLOPTV2/Release/xstelloptv2 input.QAS > log.stellopt



