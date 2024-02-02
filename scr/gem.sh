#!/bin/bash

## header  --------------------------------------------------------------------

### the one that calculates the mappability values for calling CNVs
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=22
ref_name="${1}"
read_len="${2}"

### output folder
out_dir="${base_dir}/gem/mpt"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

### input
ref_dir="${base_dir}/ref"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that calculates \
the mappability values for calling CNVs..."

cd "${out_dir}"

### make a temporary copy of the reference
ref_fa=$(find "${ref_dir}" -name "${ref_name}-genome.fa")
cp "${ref_fa}" "${out_dir}"

### gem indexer
gem-indexer -T "${n_threads}" -c dna -i "${ref_name}-genome.fa" \
-o "${ref_name}"

### calculate mappability
gem-mappability -T "${n_threads}" -I "${ref_name}.gem" \
-l "${read_len}" -o "${ref_name}"

### remove the reference copy
rm -f "${out_dir}/${ref_name}-genome.fa"


