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

db_urls=(
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20241228.tar.gz" # new standard
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20241228.tar.gz"
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20241228.tar.gz"
"https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20241228.tar.gz" # new minusB
"https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20241228.tar.gz" # new viral
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240904.tar.gz" # old standards
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz"
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240112.tar.gz"
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20231009.tar.gz"
"https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz"
)

define_vars(){
  root="/nfsd/bcb/bcbg/rossigno/metagenomics"
  ds_dir="/nfsd/bcb/bcbg/cavattoni/Algobio2023"

  PATH="/nfsd/bcb/bcbg/rossigno/tools:$PATH"

  db_archive=$(basename $db_url)

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
  results="$root/results/$db"
}

download_database(){
    echo
    echo "### PHASE 0: DOWNLOAD $db ###"
    echo

    srun wget -q -nc "$db_url"
    srun tar -xvf "$db_archive" --one-top-level
    rm -f "$db_archive"

    echo "--- Done ---"
}

run_kraken(){
  echo
  echo "### PHASE 1: RUNNING KRAKEN2 for $db ###"
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

  echo "--- Done ---"
}

evaluate(){
  ### EVALUATION ###
  echo
  echo "### PHASE 2: RUNNING EVALUATION for $db ###"
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

  srun 'NR>1{print $11}' <<< "$(evaluate_calls $nodes genus $Rground $Rground)"
  srun 'NR>1{print $11}' <<< "$(evaluate_calls $nodes genus $Sground $Sground)"
  srun 'NR>1{print $11}' <<< "$(evaluate_calls $nodes genus $Mground $Mground)"

  echo
  echo "Number of species per dataset:"
  echo

  srun 'NR>1{print $11}' <<< "$(evaluate_calls $nodes species $Rground $Rground)"
  srun 'NR>1{print $11}' <<< "$(evaluate_calls $nodes species $Sground $Sground)"
  srun 'NR>1{print $11}' <<< "$(evaluate_calls $nodes species $Mground $Mground)"

  echo
  echo "Number of genus per database ($db):"
  echo

  srun grep -c -P "\tG\t" $root/databases/$db/inspect.txt

  echo
  echo "Number of species per database ($db)"
  echo

  srun grep -c -P "\tS\t" $root/databases/$db/inspect.txt

  echo "--- Done ---"
}

generate_report(){
  echo
  echo "### PHASE 3: GENERATE KRAKEN2 REPORT for $db ###"
  echo

  srun kraken2 --report-zero-counts --db $root/databases/$db --threads 16 --paired $R1 $R2 --report $results/REPORT_REAL_$db

  srun kraken2 --report-zero-counts --db $root/databases/$db --threads 16 --paired $S1 $S2 --report $results/REPORT_SIMULATED_$db

  srun kraken2 --report-zero-counts --db $root/databases/$db --threads 16 $M --report $results/REPORT_MARINE_$db

  echo "--- Done ---"
}

echo "### JOB STARTED ###"

for db_url in "${db_urls[@]}"; do
  echo "================================================================================"

  define_vars

  cd $root/databases

  download_database

  mkdir -p "$results" && cd "$results"

  run_kraken

  evaluate

  generate_report

  echo "================================================================================"
done

echo "### JOB TERMINATED ###"
