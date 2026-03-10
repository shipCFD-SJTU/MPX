#!/bin/bash
#SBATCH --job-name=ysjsph
#SBATCH --partition=cpu
#SBATCH -n 200
#SBATCH --ntasks-per-node=40
#SBATCH --mail-type=end
#SBATCH --mail-user=13667320831@sjtu.edu.cn
#SBATCH --output=%j.out
#SBATCH --error=%j.err

ulimit -s unlimited
ulimit -l unlimited

module load openmpi/4.1.5-gcc-8.5.0

module load gsl/2.7.1-intel-2021.4.0

source ~/OpenFOAM/OpenFOAM-v2106/etc/bashrc

workhome=`pwd`
cd $workhome

decomposePar

mpirun -np 200 stratifiedFoam -parallel

reconstructPar


