#!/bin/bash

## header  --------------------------------------------------------------------

### the one that makes the fai indexes of all the genomes in the "rep" folder
### (thay are needed to build a valid vcf file from nucmer output)

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")

### input
repo_dir="${base_dir}/rep"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that \
makes the fai indexes of all the genomes in the rep folder..."

cd "${repo_dir}"
gen_fa=$(find "${repo_dir}" -name "*fa")

for ind_g in ${gen_fa}; do
  samtools faidx "${ind_g}"
done
