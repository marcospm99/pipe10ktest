#!/bin/bash
# Job Name and Files (also --job-name)
#SBATCH -J test
#Output and error (also --output, --error):
#SBATCH -o /dss/dsshome1/00/di25hek2/output/%x.%j.out
#SBATCH -e /dss/dsshome1/00/di25hek2/output/%x.%j.err
#Initial working directory (also --chdir):
#SBATCH -D /dss/dsshome1/00/di25hek2/test
#Notification and type
#SBATCH --mail-type=END
#SBATCH --mail-user=serhocal@mot.upv.es
# Wall clock limit:
#SBATCH --time=00:29:55
##SBATCH --time=47:40:55
#SBATCH --no-requeue
#Setup of execution environment
#SBATCH --export=NONE
#SBATCH --get-user-env
#SBATCH --account=pn49ya

#SBATCH --partition=test
##SBATCH --partition=fat
#Number of nodes and MPI tasks per node: 
# WARNING JUST NODES
#SBATCH --nodes=11

#####SBATCH --partition=test
####Number of nodes and MPI tasks per node: 
#### WARNING JUST NODES
####SBATCH --nodes=11
###

#SBATCH --ntasks-per-node=24

##SBATCH --ntasks-per-core=1
##SBATCH --cpus-per-task=2

##module load hdf5/1.10.2-intel-impi
### These two modules are needd by H5DF

module load slurm_setup
module load hdf5/1.10.7-intel21-impi

_NR=16
_NS=8

sed s/#_NR#/$_NR/ < parallel.tmp > tmp
sed s/#_NS#/$_NS/ < tmp > parallel.h
rm tmp
export OMP_NUM_THREADS=1

# Prepare hre.dat

OF=$(ls -1  /dss/dsshome1/00/di25hek2/work/pipe/pipe* | tail -1 | cut -f2 -d "." )
: $((OF++))

echo $OF > hre.dat
echo /dss/dsshome1/00/di25hek2/work/pipe/ >> hre.dat
echo pipe180 >> hre.dat
echo /dss/dsshome1/00/di25hek2/work/pipe/sta/pipe180_${_NR}x${_NS} >> hre.dat
procs=$((_NR*_NS))

mpiexec -np $procs ./pipe10k 
echo submitting
sbatch cpa.sh
