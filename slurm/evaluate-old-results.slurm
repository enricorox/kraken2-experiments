#!/bin/bash

#SBATCH --job-name kraken2-old
#SBATCH --output logs/kraken2-old-res.log
#SBATCH --error logs/kraken2-old-res.log
#SBATCH --mail-user enrico.rossignolo@studenti.unipd.it
#SBATCH --mail-type ALL
#SBATCH --time 14-00:00:00
#SBATCH --ntasks 1
#SBATCH --cpus-per-task=16
#SBATCH --partition allgroups
#SBATCH --mem 200G

set -e
set -x
shopt -s nullglob

define_vars(){
  root="/nfsd/bcb/bcbg/rossigno/metagenomics"
  ds_dir="/nfsd/bcb/bcbg/cavattoni/Algobio2023"

  PATH="/nfsd/bcb/bcbg/rossigno/tools:$PATH"

  db_archive=$(basename "$db_url")

  db=${db_archive%.tar.gz}
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
  results="/nfsd/bcb/bcbg/lazzaro/result"
}

evaluate_real(){
  ### EVALUATION ###
  echo
  echo "### PHASE 2: RUNNING EVALUATION ###"
  echo

  echo
  echo "Evaluation: Genus"
  echo

  srun evaluate_calls $nodes genus $Rground "$res"

  echo
  echo "Evaluation: Species"
  echo

  srun evaluate_calls $nodes species $Rground "$res"

  echo "### DONE ###"
}

evaluate_simulated(){
  ### EVALUATION ###
  echo
  echo "### PHASE 2: RUNNING EVALUATION ###"
  echo

  echo
  echo "Evaluation: Genus"
  echo

  srun evaluate_calls $nodes genus $Sground "$res"

  echo
  echo "Evaluation: Species"
  echo

  srun evaluate_calls $nodes species $Sground "$res"

  echo "### DONE ###"
}

evaluate_marine(){
  ### EVALUATION ###
  echo
  echo "### PHASE 2: RUNNING EVALUATION"
  echo

  echo
  echo "Evaluation: Genus"
  echo

  srun evaluate_calls $nodes genus $Mground "$res"

  echo
  echo "Evaluation: Species"
  echo

  srun evaluate_calls $nodes species $Mground "$res"

  echo "### DONE ###"
}

echo "### JOB STARTED ###"

define_vars

for res in "$results"/REAL*.res; do
  echo "================================================================================"

  echo "Evaluation: $res"

  evaluate_real

  #generate_report

  echo "================================================================================"
done

for res in "$results"/SMLTD*.res; do
  echo "================================================================================"

  echo "Evaluation: $res"

  evaluate_simulated

  #generate_report

  echo "================================================================================"
done

for res in "$results"/MARINE*.res; do
  echo "================================================================================"

  echo "Evaluation: $res"

  evaluate_marine

  #generate_report

  echo "================================================================================"
done

echo "### JOB TERMINATED ###"
