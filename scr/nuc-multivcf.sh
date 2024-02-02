#!/bin/bash

## header  --------------------------------------------------------------------

### the one that makes the multisample vcf file for nucmer data

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
dir_base=$(dirname "${full_dir}")
n_threads=4

### input and output folders
dir_input="${dir_base}/gen/os-par-vcf"
dir_out="${dir_base}/gen/ms-par-vcf"
if [[ -d "${dir_out}" ]]; then
  rm -rf "${dir_out}"
fi
mkdir -p "${dir_out}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that makes the multisample vcf file for nucmer data..."

cd "${dir_input}"
rm -f *tbi

### compress with bgzip all the single-sample vcf files and index with tabix
for vcf_file in $(find . -name "*vcf"); do
  bgzip -c -f "${vcf_file}" > "${vcf_file}.gz"
  tabix "${vcf_file}.gz"
done

### merge to one multisample vcf with bcftools
all_files=( $(find . -name "*vcf.gz") )
bcftools merge --output-type z --threads "${n_threads}" ${all_files[@]} \
> "${dir_out}/multis-snps.vcf.gz"

### prepare header for the next chunk
bgzip -b -c "${dir_out}/multis-snps.vcf.gz" | grep "^#" \
> "${dir_out}/multis-snps-genfix.vcf"
