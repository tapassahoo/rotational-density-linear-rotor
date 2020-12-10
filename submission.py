#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.log
#SBATCH --account=rrg-pnroy
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
python rot-dens.py 5. 20 20.561 1500 spinless 41 41
