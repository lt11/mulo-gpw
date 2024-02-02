#!/bin/bash

## header  --------------------------------------------------------------------

### the one that checks the quality of the fastq files
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=22

### input
fastq_ext=".fq.gz"
fq_dir="${base_dir}/exp"

### output
out_dir="${base_dir}/fqc"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that checks the quality of the fastq files..."

cd "${fq_dir}"

pll_check=$((pll_runs + 1))
for ind_fq in $(\ls *"${fastq_ext}"); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi
  
  fastqc "${ind_fq}" &> /dev/null &
done

wait

mv *"_fastqc"* "${out_dir}"
