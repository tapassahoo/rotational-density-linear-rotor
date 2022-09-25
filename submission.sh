#!/bin/bash
#SBATCH --job-name=test
#SBATCH --output=test.log
#SBATCH --account=rrg-pnroy
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=1GB
#SBATCH --cpus-per-task=1
export OMP_NUM_THREADS=1
python  rotational_energy_propagator.py --temperature 5.0 --trotter_number 20 --rotational_constant 20.561 --size 1500
