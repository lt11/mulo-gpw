#!/bin/bash

## header  --------------------------------------------------------------------

### the one that aligns the genomes of the parents

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
ref_name="${1}"

### output folder
out_dir="${base_dir}/gen/aln-par"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

### input
repo_dir="${base_dir}/rep"
ref_fa=$(find "${repo_dir}" -name "${ref_name}*fa")

## clmnt  ---------------------------------------------------------------------

echo "Running the one that aligns the genomes of the parents..."

cd "${out_dir}"

## run nucmer -----------------------------------------------------------------

### get all chromosomes in the reference
all_chr=$(grep "^>" "${ref_fa}" | sed 's|^.||')
for ind_c in ${all_chr}; do
  ref_chr="${ref_name}-${ind_c}.fa"
  samtools faidx "${ref_fa}" "${ind_c}" > "${ref_chr}"
done

### get the patchwork genomes
patchwork_genomes=$(find "${repo_dir}" -name "*fa" | grep -v "${ref_name}")

for ind_g in ${patchwork_genomes}; do
  ### this will always provide e.g. Y12 as prefix whatever the name of the
  ### fasta, e.g. Y12-hc-genome.fa (the current, PanSN-spec aware, standard)
  ### or Y12-genome.fa (the old standard)
  pref_pword=$(basename "${ind_g}" | cut -d "-" -f 1)
  pref_align="${ref_name}-${pref_pword}"
  ### cleaning output files
  if [[ -f "${pref_align}-var.txt" ]]; then
    rm "${pref_align}-var.txt"
  fi
  if [[ -f "${pref_align}-coords.txt" ]]; then
    rm "${pref_align}-coords.txt"
  fi
  for ind_c in ${all_chr}; do
    ref_chr="${ref_name}-${ind_c}.fa"
    pwork_chr="${pref_pword}-${ind_c}.fa"
    samtools faidx "${ind_g}" "${ind_c}" > "${pwork_chr}"
    nucmer --prefix "${pref_align}" \
    "${ref_chr}" "${pwork_chr}"
    show-snps -CrTH "${pref_align}.delta" \
    >> "${pref_align}-var.txt"
    show-coords -TH "${pref_align}.delta" \
    >> "${pref_align}-coords.txt"
  done
done

### remove the fastas of the chromosomes
rm -f *chr*.fa
