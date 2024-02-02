#!/bin/bash

## header  --------------------------------------------------------------------

### the one which intersects the variants from bcftools and freebayes
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname $0) && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=22
pll_samples=2

### output folder
out_dir="${base_dir}/var-calls/int-markers"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir "${out_dir}"

### input folder
var_dir="${base_dir}/var-calls/raw-markers"

## clmnt  ---------------------------------------------------------------------

echo "Running the one which intersect the variants \
from bcftools and freebayes..."

cd "${var_dir}"
for ind_bcft_file in $(ls *"-bcft.vcf.gz"); do
  if (( cnt_p % pll_samples == 0 )); then
    wait
  fi
  ((cnt_p++))
  (
  ind_freeb_file=$(echo "${ind_bcft_file}" | sed 's|bcft|freeb|')
  ind_int_file=$(echo "${ind_bcft_file}" | sed 's|bcft|int|')
  vcf-isec -n +2 "${ind_bcft_file}" "${ind_freeb_file}" | bgzip \
  > "${out_dir}/${ind_int_file}" 
  tabix -f -p vcf "${out_dir}/${ind_int_file}"
  ) &
done

wait


