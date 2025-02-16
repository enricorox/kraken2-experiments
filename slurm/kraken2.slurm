#!/bin/bash

#SBATCH --job-name kraken2
#SBATCH --output logs/kraken2.log
#SBATCH --error logs/kraken2.log
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH --time 1-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=16
#SBATCH --partition allgroups
#SBATCH --mem 200G

set -e

root="/nfsd/bcb/bcbg/rossigno/metagenomics"
ds_dir="/nfsd/bcb/bcbg/cavattoni/Algobio2023"

PATH="/nfsd/bcb/bcbg/rossigno/tools:$PATH"

db="Standard20240605" # change database here
nodes="/nfsd/bcb/bcbg/rossigno/metagenomics/reads/nodes.dmp"

# real dataset and ground truth
R1="$ds_dir/SRR1804065/SRR1804065_1.filtr.fq"
R2="$ds_dir/SRR1804065/SRR1804065_2.filtr.fq"
Rground="/nfsd/bcb/bcbg/cavattoni/Algobio2023/SRR1804065/SRR1804065_truth.tsv"

# Simulated dataset and ground truth
S1="$ds_dir/Sim/all_250000_1.fq.gz"
S2="$ds_dir/Sim/all_250000_2.fq.gz"
Sground="/nfsd/bcb/bcbg/cavattoni/Algobio2023/Sim/truth_all_250000.tsv"

# Marine dataset and ground truth
M="$ds_dir/Marine/marmgCAMI2_long_read_sample_0_reads.fq"
Mground="/nfsd/bcb/bcbg/cavattoni/Algobio2023/Marine/marmgCAMI2_long_read_sample_0_reads.fa.blastn.txt"

# Results directory
results="$root/results/$db"

mkdir -p $results && cd $results

### KRAKEN2 MAIN ###

cd $results

run_kraken(){
  echo
  echo "### PHASE 1: RUNNING KRAKEN ###"
  echo

  #REAL DATASET
  echo
  echo "Running Kraken2 with $db on real dataset"
  echo

  srun kraken2 --db $root/databases/$db --threads 16 --paired $R1 $R2 > $results/REAL_$db.out

  #SIMULATED DATASET
  echo
  echo "Running Kraken2 with $db on simulated dataset"
  echo

  srun kraken2 --db $root/databases/$db --threads 16 --paired $S1 $S2 > $results/SMLTD_$db.out

  #MARINE DATASET
  echo
  echo "Running Kraken2 with $db on marine dataset"
  echo

  srun kraken2 --db $root/databases/$db --threads 16 $M > $results/MARINE_$db.out

  # refining results
  srun cut -f 2-3 $results/REAL_$db.out > $results/REAL_$db.res
  srun cut -f 2-3 $results/SMLTD_$db.out > $results/SMLTD_$db.res
  srun cut -f 2-3 $results/MARINE_$db.out > $results/MARINE_$db.res
}

evaluate(){
  ### EVALUATION ###
  echo
  echo "### PHASE 2: RUNNING EVALUATION ###"
  echo

  echo
  echo "Evaluation: Genus of $db"
  echo

  srun evaluate_calls $nodes genus $Rground $results/REAL_$db.res
  srun evaluate_calls $nodes genus $Sground $results/SMLTD_$db.res
  srun evaluate_calls $nodes genus $Mground $results/MARINE_$db.res

  echo
  echo "Evaluation: Species of $db"
  echo

  srun evaluate_calls $nodes species $Rground $results/REAL_$db.res
  srun evaluate_calls $nodes species $Sground $results/SMLTD_$db.res
  srun evaluate_calls $nodes species $Mground $results/MARINE_$db.res

  echo
  echo "Number of genus per dataset (see DistinctTPAtRank)"
  echo

  srun evaluate_calls $nodes genus $Rground $Rground
  srun evaluate_calls $nodes genus $Sground $Sground
  srun evaluate_calls $nodes genus $Mground $Mground

  echo
  echo "Number of species per dataset (see DistinctTPAtRank)"
  echo

  srun evaluate_calls $nodes species $Rground $Rground
  srun evaluate_calls $nodes species $Sground $Sground
  srun evaluate_calls $nodes species $Mground $Mground

  echo
  echo "Number of genus per database ($db)"
  echo

  srun grep -c -P "\tG\t" $root/databases/$db/inspect.txt

  echo
  echo "Number of species per database ($db)"
  echo

  srun grep -c -P "\tS\t" $root/databases/$db/inspect.txt
}

echo "### JOB STARTED ###"

#run_kraken

evaluate

echo "### JOB TERMINATED ###"
