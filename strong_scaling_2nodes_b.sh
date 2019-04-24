#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=cardiac-sim
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=16
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --partition=short
##SBATCH --exclusive
##SBATCH --constraint=e52695v4,32cpu
##SBATCH --time=30:00
#SBATCH --output=cardiacsim-2nodes.out
##SBATCH --mail-type=ALL

################################################################################
################################################################################

## Load openmpi version 3.0.0
echo "Loading openmpi module ..."
module load mpich/3.2.1

## Load GCC-7.2.1
echo "Loading GCC module ..."
module load gcc/7.3.0

echo ""
echo "======================================================================================"
env
echo "======================================================================================"
echo ""

# Set stack size to unlimited
echo "Setting stack size to unlimited..."
ulimit -s unlimited
ulimit -l unlimited
ulimit -a
echo

cd parallel/A2/cardiacsim

echo "With 1D geometry :"
mpirun -np 32 ./cardiacsim_only_mpi -x 1 -y 32 -n 1024 -t 100

echo "With 2D geometry :"
mpirun -np 32 ./cardiacsim_only_mpi -x 4 -y 8 -n 1024 -t 100

#....
echo "Finished with execution!"
