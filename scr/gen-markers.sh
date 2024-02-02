#!/bin/bash

## header  --------------------------------------------------------------------

### the one that genotypes the samples at markers positions

## settings  ------------------------------------------------------------------

### dev 
# full_dir="/home/ltattini/prog/evolrec-ae/gpw-pm-1/scr"
# base_dir="/home/ltattini/prog/evolrec-ae/gpw-pm-1"
# n_ploidy="2"
full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
n_ploidy="${1}"
exp_des="${2}"
n_threads=8
pll_runs=4

### input
map_dir="${base_dir}/map-sr"
flt_dir="${base_dir}/var-calls/flt-markers"
ref_file=$(find "${base_dir}/ref" -name "*-genome.fa")

### output
out_dir="${base_dir}/var-calls/gen-markers"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that genotypes the samples at markers positions..."

### compress marker files and retrieve the ID of the parent strain
cd "${flt_dir}"
pll_check=$((pll_runs + 1))
for vcf_file in $(find . -name "*vcf"); do
  bgzip -c -f "${vcf_file}" > "${vcf_file}.gz"
  tabix "${vcf_file}.gz"
  mrk_id=$(echo "${vcf_file}" | cut -d "-" -f 2 | cut -d "." -f 1)
  for ind_bam in $(find "${map_dir}" -name "*bam"); do
    ### parallel samples
    ((cnt_p++))
    if (( cnt_p % pll_check == 0 )); then
      wait -n
      cnt_p=$(( pll_check - 1 ))
    fi
    
    (
    ### genotyping
    samp_id=$(basename "${ind_bam}" | cut -d "-" -f 1)
    out_file="${out_dir}/${samp_id}-${mrk_id}.vcf.gz"
    bcftools mpileup --threads "${n_threads}" -T "${vcf_file}.gz" \
    -min-MQ5 -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP \
    --skip-indels -D -f "${ref_file}" "${ind_bam}" | \
    bcftools call --threads "${n_threads}" --ploidy "${n_ploidy}" -c -O z \
    > "${out_file}"
    tabix -f -p vcf "${out_file}"
    ) &
  done
done

wait

### obsolete since now we detect markers from the assemblies:
### get the experimental design, if it is with not plus/minus samples
### (hence it is an evolved/ancestor desing)
### we need a copy of the ancestor vcf file in the gen-markers folder
# exp_des=$(grep "^exp_design=" mulo-*sh | cut -d "=" -f 2 | sed 's|"||g')
# if [[ "${exp_des}" == "BH" ]]; then
#   mrk_id=$(basename "${mrk_file}" | cut -d "-" -f 1)
#   cp "${mrk_file}" "${out_dir}/${mrk_id}.vcf.gz"
#   tabix -f -p vcf "${out_dir}/${mrk_id}.vcf.gz"
# fi

# the definition in the header:
### the experimental design ("exp_design" variable) can be "A"
### (for plus/minus phenotyped samples)
### or "BH" (for evolved/ancestral samples)
# exp_design="A"
