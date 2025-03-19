#!/bin/bash

inspects=(
https://genome-idx.s3.amazonaws.com/kraken/standard_20221209/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_20221209/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20221209/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_20220926/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_20220908/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20220908/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_20220607/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_20220607/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20220607/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_20210517/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_20210517/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_20201202/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_20210127/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20210127/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_20200919/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_20200919/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspfp_20200919/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/viral_20220908/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/minusb_20220926/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_08gb_20220926/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_16gb_20220926/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_08gb_20220908/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/pluspf_16gb_20220908/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_20230605/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_08gb_20230605/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_16gb_20230605/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_20230314/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_08gb_20230314/inspect.txt
https://genome-idx.s3.amazonaws.com/kraken/standard_16gb_20230314/inspect.txt
)

count(){
  echo
    echo "Number of genus per database ($db):"
    echo

    grep -c -P "\tG\t" "inspect.txt"

    echo
    echo "Number of species per database ($db)"
    echo

    grep -c -P "\tS\t" "inspect.txt"

    echo "### DONE ###"
}

echo "### JOB STARTED ###"

for inspect in "${inspects[@]}"; do
  db=$(basename "${inspect%/inspect.txt}")
  echo "Downloading inspect.txt for $db"
  wget -q "$inspect"
  count
  rm -f inspect.txt
done

echo "### JOB TERMINATED ###"
