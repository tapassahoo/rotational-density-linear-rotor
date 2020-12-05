#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.log
#SBATCH --account=rrg-pnroy
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
python rot-dens.py 10. 100 60.853 1500 spinless 51 51
