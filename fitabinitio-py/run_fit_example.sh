#!/bin/bash
 
# set the account to be used for the job
#SBATCH --account=free
 
# set name of job
#SBATCH --job-name=FT
#SBATCH --output=FT.%j.o
#SBATCH --error=FT.%j.e
 
# set the number of nodes and partition
#SBATCH --nodes=4
#####BATCH --ntasks-per-node=4
####SBATCH --partition=batch-64gb
#SBATCH --partition=batch-devel
#SBATCH --qos=devel

 
# set max wallclock time
#SBATCH --time=00:15:00
 
 
# Load dependant modules
export PYTHONPATH=/home/c/cor22/scratch/ppfit_full_mpi/
export I_MPI_FABRICS=dapl,ofa,tcp,tmi,ofi
export I_MPI_FALLBACK=yes

module load anaconda3/2.5.0
module load intel/mpi/64/5.1.3.210

# automatically set up host file:

srun hostname -s | sort -u > slurm.hosts
echo $SLURM_NNODES

for ((i=0; i<$SLURM_NNODES; i++)); do
  echo ":16">>proc
done
paste slurm.hosts proc > temp
awk '{print $1$2}' temp > slurm.hosts


rm temp
rm proc

#mpi4py:
mpiexec.hydra -genv I_MPI_FABRICS=dapl,ofa,tcp,tmi,ofi --rr --map-by node --bind-to node -f slurm.hosts -np 33 python3 ./fitabinitio.py
#serial & pool:
#./fitabinitio.py



