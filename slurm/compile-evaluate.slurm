#!/bin/bash

#SBATCH --job-name compile-eval
#SBATCH --output compile-evaluate.log
#SBATCH --error compile-evaluate.log
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH --time 1-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=1
#SBATCH --partition allgroups
#SBATCH --mem 4G

set -e

cd "/nfsd/bcb/bcbg/rossigno/metagenomics/kraken2-experiments"

srun cmake .

srun make

ln -s "$(realpath evaluate_calls)" "/nfsd/bcb/bcbg/rossigno/tools/evaluate_calls"
