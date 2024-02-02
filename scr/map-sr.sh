#!/bin/bash

## header  --------------------------------------------------------------------

### the one that maps short-reads and fixes the bam file (fixmate and markdup)

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=22
pll_runs=4
ref_name="${1}"

### output folder
out_dir="${base_dir}/map-sr"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that maps short-reads and \
fixes the bam file (fixmate and markdup)..."

cd "${base_dir}"
ref_path=$(find "${base_dir}/ref" -name "${ref_name}*fa")

all_fqs=$(ls "exp/"*"gz" | cut -d "-" -f 1 | cut -d "/" -f 2 | sort | uniq)
pll_check=$((pll_runs + 1))
for ind_e in ${all_fqs}; do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi
  
  (
  ### mapping
  bwa mem -M -t "${n_threads}" "${ref_path}" \
  "exp/${ind_e}-R1.fq.gz" \
  "exp/${ind_e}-R2.fq.gz" \
  > "${out_dir}/${ind_e}-${ref_name}.sam"

  ### fix mates
  samtools fixmate -O bam,level=1 -@ "${n_threads}" \
  -m "${out_dir}/${ind_e}-${ref_name}.sam" \
  "${out_dir}/${ind_e}-${ref_name}.bam"

  ### sorting
  samtools sort -O bam,level=1 -@ "${n_threads}" \
  -o "${out_dir}/${ind_e}-${ref_name}-srt.bam" \
  "${out_dir}/${ind_e}-${ref_name}.bam"

  ### marking duplicates
  samtools markdup -O bam,level=1 -@ "${n_threads}" \
  "${out_dir}/${ind_e}-${ref_name}-srt.bam" \
  "${out_dir}/${ind_e}-${ref_name}-srt-mdp.bam"

  ### indexing
  samtools index -@ "${n_threads}" \
  "${out_dir}/${ind_e}-${ref_name}-srt-mdp.bam"

  ### cleaning
  rm -f "${out_dir}/${ind_e}-${ref_name}.sam"
  rm -f "${out_dir}/${ind_e}-${ref_name}.bam"
  rm -f "${out_dir}/${ind_e}-${ref_name}-srt.bam"
  ) &
done

wait
