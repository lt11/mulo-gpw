#!/bin/bash

## header  --------------------------------------------------------------------

### the one that calls markers using bcftools and freebayes
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname $0) && pwd)
base_dir=$(dirname "${full_dir}")
n_ploidy="${1}"
all_ancestors="${2}"
n_threads=22
pll_samples=2

### input
data_ext="bam"
map_dir="${base_dir}/map-sr"

### output folder
out_dir="${base_dir}/var-calls/raw-markers"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that calls markers using bcftools and freebayes..."

cd "${map_dir}"

for ind_samp in ${all_ancestors}; do
  ind_map=$(ls "${ind_samp}-"*"${data_ext}")
  ### parallel samples
  if (( cnt_p % pll_samples == 0 )); then
    wait
  fi
  ((cnt_p++))
  (
  ref_path=$(samtools view -H "${ind_map}" | grep "^@PG" \
  | head -1 | cut -f 5 | cut -d " " -f 6)
  ref_name=$(basename "${ref_path}" | cut -d "-" -f 1)
  samp_name=$(echo "${ind_map}" | cut -d "-" -f 1)
  
  ### bcftools: the output is compressed
  (
  out_bcft="${out_dir}/${samp_name}-${ref_name}-bcft.vcf.gz"
  bcftools mpileup --threads "${n_threads}" \
  --max-depth 2000 -Ou -f "${ref_path}" "${ind_map}" \
  | bcftools call --threads "${n_threads}" \
  -mv -O z --ploidy "${n_ploidy}" > "${out_bcft}"
  tabix -f -p vcf "${out_bcft}"
  ) &
  
  ### freebayes: the output is NOT compressed
  out_freeb="${out_dir}/${samp_name}-${ref_name}-freeb.vcf"
  freebayes -p "${n_ploidy}" --fasta-reference "${ref_path}" "${ind_map}" \
  --limit-coverage 2000 > "${out_freeb}"
  ### fix sample name
  sed "s|unknown|${ind_map}|" "${out_freeb}" > "${out_freeb}.tmp"
  \mv -f "${out_freeb}.tmp" "${out_freeb}"
  bgzip "${out_freeb}"
  tabix -f -p vcf "${out_freeb}.gz"
  ) &
done

wait


