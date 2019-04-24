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
#SBATCH --ntasks-per-node=32
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

echo "Serial version ..."
./cardiacsim_serial -n 1024 -t 100

echo "With 1D geometry :"
echo "1 MPI"
mpirun -np 1 ./cardiacsim_only_mpi -x 1 -y 1 -n 1024 -t 100
echo "2 MPI"
mpirun -np 2 ./cardiacsim_only_mpi -x 1 -y 2 -n 1024 -t 100
echo "4 MPI"
mpirun -np 4 ./cardiacsim_only_mpi -x 1 -y 4 -n 1024 -t 100
echo "8 MPI"
mpirun -np 8 ./cardiacsim_only_mpi -x 1 -y 8 -n 1024 -t 100
echo "16 MPI"
mpirun -np 16 ./cardiacsim_only_mpi -x 1 -y 16 -n 1024 -t 100
echo "32 MPI"
mpirun -np 32 ./cardiacsim_only_mpi -x 1 -y 32 -n 1024 -t 100

echo "With 2D geometry :"
echo "1 MPI"
mpirun -np 1 ./cardiacsim_only_mpi -x 1 -y 1 -n 1024 -t 100
echo "2 MPI"
mpirun -np 2 ./cardiacsim_only_mpi -x 1 -y 2 -n 1024 -t 100
echo "4 MPI"
mpirun -np 4 ./cardiacsim_only_mpi -x 2 -y 2 -n 1024 -t 100
echo "8 MPI"
mpirun -np 8 ./cardiacsim_only_mpi -x 2 -y 4 -n 1024 -t 100
echo "16 MPI"
mpirun -np 16 ./cardiacsim_only_mpi -x 4 -y 4 -n 1024 -t 100
echo "32 MPI"
mpirun -np 32 ./cardiacsim_only_mpi -x 4 -y 8 -n 1024 -t 100

#....
echo "Finished with execution!"
