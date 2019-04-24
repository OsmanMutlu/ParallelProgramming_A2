#!/bin/bash
#
# You should only work under the /scratch/users/<username> directory.
#
# Example job submission script
#
# -= Resources =-
#
#SBATCH --job-name=cardiac-sim
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=users
#SBATCH --qos=users
#SBATCH --partition=short
##SBATCH --exclusive
##SBATCH --constraint=e52695v4,32cpu
##SBATCH --time=30:00
#SBATCH --output=cardiacsim-%j.out
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
echo "N = 1024:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 1024 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 1024 -t 100 -k
echo "---------------------------------------------------"
echo "N = 768:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 768 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 768 -t 100 -k
echo "---------------------------------------------------"
echo "N = 512:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 512 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 512 -t 100 -k
echo "---------------------------------------------------"
echo "N = 384:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 384 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 384 -t 100 -k
echo "---------------------------------------------------"
echo "N = 256:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 256 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 256 -t 100 -k

echo "====================================================="

echo "With 2D geometry :"
echo "N = 1024:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 1024 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 1024 -t 100 -k
echo "---------------------------------------------------"
echo "N = 768:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 768 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 768 -t 100 -k
echo "---------------------------------------------------"
echo "N = 512:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 512 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 512 -t 100 -k
echo "---------------------------------------------------"
echo "N = 384:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 384 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 384 -t 100 -k
echo "---------------------------------------------------"
echo "N = 256:"
echo "With Comm:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 256 -t 100
echo "With Comm Disabled:"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 256 -t 100 -k

#....
echo "Finished with execution!"
