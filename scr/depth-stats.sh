#!/bin/bash

## header  --------------------------------------------------------------------

### the one that makes the statistics of the coverage
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=22

### output folder
out_dir="${base_dir}/cov"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

## clmnt  ---------------------------------------------------------------------

cd "${base_dir}"

pll_check=$((pll_runs + 1))
for ind_m in $(ls "map-sr/"*"bam"); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi
  
  (
  map_name=$(echo "${ind_m}" | cut -d "/" -f 2 | cut -d "-" -f 1,2)
  samtools depth "${ind_m}" > "${out_dir}/${map_name}-depth.txt"
  ref_len=$(samtools view -H "${ind_m}" | \
  grep "^@SQ" | cut -f 3 | cut -d ":" -f 2 | \
  awk '{sum += $1} END {print sum}')
  awk -v RL="${ref_len}" 'BEGIN {FS="\t"} {sum+=$3; sumsq+=$3*$3} \
  END {print "Covered sites=" NR; print "Reference covered (%)=" 100*NR/RL; \
  print "Reference size=" RL; print "Mean (depth)=",sum/NR; \
  print "SD (depth)=",sqrt(sumsq/NR - (sum/NR)**2)}' \
  "${out_dir}/${map_name}-depth.txt" > "${out_dir}/${map_name}-mean.txt"
  ) &
done

wait

rm -f "cov/"*"depth.txt"



