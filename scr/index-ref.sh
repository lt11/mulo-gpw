#!/bin/bash

## header  --------------------------------------------------------------------

### the one that creates the ref folder

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
rep_dir="${base_dir}/rep"
ref_name="${1}"

### output folder
out_dir="${base_dir}/ref"
if [[ ! -d "${out_dir}" ]]; then mkdir -p "${out_dir}"; fi

## clmnt  ---------------------------------------------------------------------

echo "Running the one that creates the ref folder..."

### copy the reference and make indexes for BWA
cp "${rep_dir}/${ref_name}-genome.fa" "${out_dir}"
ref_path=$(find "${out_dir}" -name "${ref_name}*fa")
bwa index "${ref_path}"


