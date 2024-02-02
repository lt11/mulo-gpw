#!/bin/bash

## header  --------------------------------------------------------------------

### the one that filters the markers detected

## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname $0) && pwd)
base_dir=$(dirname "${full_dir}")
n_threads=22
pll_samples=2

### input
int_dir="${base_dir}/var-calls/int-markers"

### output
out_dir="${base_dir}/var-calls/flt-markers"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir "${out_dir}"
stat_file="${out_dir}/qual-stat.txt"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that filters the markers detected..."

cd "${int_dir}"
for ind_int_file in $(find . -name *"-int.vcf.gz"); do
  if (( cnt_p % pll_samples == 0 )); then
    wait
  fi
  ind_flt_file=$(echo "${ind_int_file}" | sed 's|int|flt|' | sed 's|.gz$||')
  samp_id=$(echo "${ind_int_file}" | cut -d "-" -f 1)
  echo "Processing sample: ${samp_id}"
  ### print the header to output file
  zgrep "^#" "${ind_int_file}" \
  > "${out_dir}/${ind_flt_file}"
  echo "Number of intersected variants: \
  $(zgrep -v "^#" "${ind_int_file}" | wc -l)"
  ### filter out weird chromosomes and indels
  zgrep -v "^#" "${ind_int_file}" | \
  grep -v "mating_type_region" | \
  grep -v "chr_II_telomeric_gap" | \
  grep -v "mitochondrial" | \
  grep -v "INDEL" \
  > "${out_dir}/${samp_id}.tmp"
  echo "Markers before QUAL filtering: \
  $(wc -l ${out_dir}/${samp_id}.tmp | cut -d " " -f 1)" > "${stat_file}"
  ### calculate mean QUAL minus one (population) standard deviation
  qual_threshold=$(awk 'BEGIN {FS="\t"} \
  {x+=$6;y+=$6^2} END {print x/NR - sqrt(y/NR-(x/NR)^2)}' \
  "${out_dir}/${samp_id}.tmp")
  echo "Quality threshold: ${qual_threshold}" >> "${stat_file}"
  ### filter for QUAL, filter out multi-allelic markers (if any),
  ### and append the results
  awk -v val="${qual_threshold}" 'BEGIN {FS="\t"} {if ($6 > val \
  && length($5) == 1 \
  && $5 ~ /[A,T,C,G]/ ) \
  print $0}' "${out_dir}/${samp_id}.tmp" >> "${out_dir}/${ind_flt_file}"
  echo "High-quality markers: \
  $(grep -v "^#" ${out_dir}/${ind_flt_file} | wc -l)" >> "${stat_file}"
  ### clean temporary file
  rm "${out_dir}/${samp_id}.tmp"
  bgzip "${out_dir}/${ind_flt_file}"
done
