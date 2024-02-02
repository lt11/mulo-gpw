#!/bin/bash

## header  --------------------------------------------------------------------

### the one for the allele frequency plots (smooth-plot.R) 
### and the statistical model (BRM.R)

## settings  ------------------------------------------------------------------

### arguments
plus_samp="${1}"
minu_samp="${2}"
exp_des="${3}"
### working folder
full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
### dev 
# full_dir="/home/ltattini/prog/evolrec-ae/gpw-pm-1/scr"
# base_dir="/home/ltattini/prog/evolrec-ae/gpw-pm-1"
# plus_samp="D1 D4 D7 D10"
# minu_samp="D2 D5 D8 D11"
# exp_des="A"

### input
gen_dir="${base_dir}/var-calls/gen-markers"
lib_dir="${base_dir}/lib"
ref_dir="${base_dir}/ref"
ref_fai=$(find "${ref_dir}" -name "*-genome.fa.fai")

### output
mrg_dir="${base_dir}/var-calls/mrg-markers"
if [[ -d "${mrg_dir}" ]]; then rm -rf "${mrg_dir}"; fi
mkdir -p "${mrg_dir}"
brm_out_dir="${base_dir}/allele-shift/brm-output"
if [[ ! "${brm_out_dir}" ]]; then rm -rf "${brm_out_dir}"; fi
mkdir -p "${brm_out_dir}"
### this is used by smooth-plot.R
mrg_multi_dir="${base_dir}/var-calls/mrg-markers-multialleles"
if [[ ! "${mrg_multi_dir}" ]]; then rm -rf "${mrg_multi_dir}"; fi
mkdir -p "${mrg_multi_dir}"
### this is used by smooth-plot.R
difalt_dir="${base_dir}/var-calls/gen-markers-diffalt"
if [[ ! "${difalt_dir}" ]]; then rm -rf "${difalt_dir}"; fi
mkdir -p "${difalt_dir}"
### this is used by smooth-plot.R
brm_in_dir="${base_dir}/allele-shift/brm-input"
if [[ -d "${brm_in_dir}" ]]; then rm -rf "${brm_in_dir}"; fi
mkdir -p "${brm_in_dir}"
### this is used by smooth-plot.R
brm_af="${base_dir}/allele-shift/af"
if [[ ! "${brm_af}" ]]; then rm -rf "${brm_af}"; fi
mkdir -p "${brm_af}"
### this is used by smooth-plot.R
afdif_dir="${base_dir}/allele-shift/plots/af-diff"
if [[ ! "${afdif_dir}" ]]; then rm -rf "${afdif_dir}"; fi
mkdir -p "${afdif_dir}"
### this is used by smooth-plot.R
af_dir="${base_dir}/allele-shift/plots/af"
if [[ ! "${af_dir}" ]]; then rm -rf "${af_dir}"; fi
mkdir -p "${af_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one for the plots (smooth-plot.R) \
and the statistical model (BRM.R)..."

### point to the correct configuration file of BRM
if [[ "${exp_des}" == "A" ]]; then
  conf_file="${lib_dir}/da-brm-config.txt"
elif [[ "${exp_des}" == "BH" ]]; then
  conf_file="${lib_dir}/db-brm-config.txt"
elif [[ "${exp_des}" == "BL" ]]; then
  conf_file="${lib_dir}/db-brm-config.txt"
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

### get all parents labels from the folder of filtered markers
### in a very funny way
dir_flt="${base_dir}/var-calls/flt-markers"
all_parents=$(find "${dir_flt}" -name "*.vcf.gz.tbi" | \
rev | cut -d "-" -f 1 | rev | cut -d "." -f 1)

### a loop on the sample/control (can be evolved/non-evolved or plus/minus)
### they are sorted in this order in the merged file
plus_dim=$(echo "${#plus_arr[@]}")
for (( ind_i=0; ind_i<plus_dim; ind_i++ )); do
  for ind_par in ${all_parents}; do
    ### output file
    mrg_file=$(echo "${plus_arr[ind_i]}-${minu_arr[ind_i]}-${ind_par}"-mrg.txt.gz)
    ### make the merge of the two to highlight missing genotypes
    bcftools merge --output-type v \
    "${gen_dir}/${plus_arr[ind_i]}-${ind_par}.vcf.gz" \
    "${gen_dir}/${minu_arr[ind_i]}-${ind_par}.vcf.gz" | \
    grep -v "^##" | bgzip > "${mrg_dir}/${mrg_file}"
    
    ### calculate AF for BRM (without mitochondrial data)
    ### and plot the AF and the smoothed values
    Rscript smooth-plot.R "${mrg_dir}/${mrg_file}"
    
    ### run BRM with the correct model
    bsa_in=$(echo "${mrg_file}" | sed 's|-mrg.txt.gz|-bsa.txt|')
    Rscript BRM.R "${conf_file}" \
    "${ref_dir}/chrom-len-brm.txt" \
    "${brm_in_dir}/${bsa_in}"
    
    blocks_file="${brm_out_dir}/${plus_arr[ind_i]}-${minu_arr[ind_i]}-${ind_par}-blocks.txt"
    qtl_file="${brm_out_dir}/${plus_arr[ind_i]}-${minu_arr[ind_i]}-${ind_par}-qtl.txt"
    mv result/result1.xls "${blocks_file}"
    mv result/result2.xls "${qtl_file}"
    rm -r result
  done
done
