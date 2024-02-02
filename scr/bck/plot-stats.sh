#!/bin/bash

## header  --------------------------------------------------------------------

### the one for the allele frequency plots and the statistical model

## settings  ------------------------------------------------------------------

### parameters
n_threads=22
pll_samples=2
### arguments
plus_samp="${1}"
minu_samp="${2}"
exp_des="${3}"
### working folder
full_dir=$(cd $(dirname $0) && pwd)
base_dir=$(dirname "${full_dir}")

### input
gen_dir="${base_dir}/var-calls/gen-markers"
aux_dir="${base_dir}/aux"
ref_dir="${base_dir}/ref"
ref_fai=$(find "${ref_dir}" -name *-genome.fa.fai)

### output
mrg_dir="${base_dir}/var-calls/mrg-markers"
if [[ -d "${mrg_dir}" ]]; then rm -rf "${mrg_dir}"; fi
mkdir "${mrg_dir}"
brm_in_dir="${base_dir}/allele-shift/brm-input"
if [[ ! -d "${brm_in_dir}" ]]; then mkdir -p "${brm_in_dir}"; fi
brm_out_dir="${base_dir}/allele-shift/brm-output"
if [[ ! -d "${brm_out_dir}" ]]; then mkdir -p "${brm_out_dir}"; fi

## clmnt  ---------------------------------------------------------------------

echo "Running the one for the plots and the statistical model..."

### point to the correct configuration file of brm
if [[ "${exp_des}" == "A" ]]; then
  conf_file="${aux_dir}/da-brm-config.txt"
elif [[ "${exp_des}" == "BH" ]]; then
  conf_file="${aux_dir}/db-brm-config.txt"
else
  echo "Invalid experimental design provided"
  echo "Accepted values: A, BH"
  echo "Provided value:" "${exp_des}"
  exit 128
fi

### make the file with chromosome sizes without the
### mitochondrial genome chrMT
grep -v "chrMT" "${ref_fai}" | cut -f 1,2 > "${ref_dir}/chrom-len-brm.txt"

### make the array of the sample IDs
read -a plus_arr <<< "${plus_samp}"
read -a minu_arr <<< "${minu_samp}"

### a loop on the sample/control (can be evolved/non-evolved or plus/minus
### they are sorted in this order in the merged file
plus_dim=$(echo "${#plus_arr[@]}")
for (( ind_i=0; ind_i<plus_dim; ind_i++ )); do
  ### output file
  out_name=$(echo "${plus_arr[ind_i]}-${minu_arr[ind_i]}"-mrg.txt.gz)
  ### make the merge of the two to highlight missing genotypes
  bcftools merge --output-type v \
  "${gen_dir}/${plus_arr[ind_i]}.vcf.gz" \
  "${gen_dir}/${minu_arr[ind_i]}.vcf.gz" | \
  grep -v "^##" | bgzip > "${mrg_dir}/${out_name}"

  ### calculate AF for BRM (without mitochondrial data)
  ### and plot the AF and the smoothed values
  Rscript smooth-plot.R "${mrg_dir}/${out_name}"

  ### run BRM with the correct model
  bsa_in=$(echo "${out_name}" | sed 's|-mrg.txt.gz|-bsa.txt|')
  Rscript BRM.R "${conf_file}" \
  "${ref_dir}/chrom-len-brm.txt" \
  "${brm_in_dir}/${bsa_in}"

  blocks_file="${brm_out_dir}/${plus_arr[ind_i]}-${minu_arr[ind_i]}-blocks.txt"
  qtl_file="${brm_out_dir}/${plus_arr[ind_i]}-${minu_arr[ind_i]}-qtl.txt"
  mv result/result1.xls "${blocks_file}"
  mv result/result2.xls "${qtl_file}"
  rm -r result
done
