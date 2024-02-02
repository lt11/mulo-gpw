## header  --------------------------------------------------------------------

### the one the calculates the CNV profiles with Control-FREEC
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
pll_runs=4

### path to template files
lib_dir="${base_dir}/lib"
tplt_file="${lib_dir}/freec-template.txt"

### input data
map_dir="${base_dir}/map-sr"
data_ext="bam"

### output folders
out_dir="${base_dir}/var-calls/cnv"
if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
mkdir -p "${out_dir}"

config_dir="${base_dir}/var-calls/config-cnv"
if [[ -d "${config_dir}" ]]; then rm -rf "${config_dir}"; fi
mkdir -p "${config_dir}"

## clmnt  ---------------------------------------------------------------------

echo "Running the one the calculates the CNV profiles with Control-FREEC..."

### retrieve ploidy from template for plotting
n_ploidy=$(grep "ploidy =" "${tplt_file}" | cut -d "=" -f 2 | sed 's|^ ||')

if [[ ! -s "${tplt_file}" ]]; then
  echo "Missing template file"
  exit 1
fi

echo "Reading ploidy from ${tplt_file}"
echo "Using ploidy = ${n_ploidy}"

### enter mappings folder
cd "${map_dir}"

pll_check=$((pll_runs + 1))
for ind_map in $(ls *"${data_ext}"); do
  ### parallel samples
  ((cnt_p++))
  if (( cnt_p % pll_check == 0 )); then
    wait -n
    cnt_p=$(( pll_check - 1 ))
  fi
  
  (
  ### retrieve reference infos
  ref_path=$(samtools view -H "${ind_map}" | grep "^@PG" \
  | head -1 | cut -f 5 | cut -d " " -f 6)
  ref_name=$(basename "${ref_path}" | cut -d "-" -f 1)
  samp_name=$(echo "${ind_map}" | cut -d "-" -f 1)
  
  ### sample-specific output
  samp_out_dir="${out_dir}/${samp_name}"
  mkdir "${samp_out_dir}"
  
  ### suffix of bam files
  bam_sfx=$(echo "${ind_map}" | cut -d "-" -f 3-)
  
  ### make the config file
  my_config="${config_dir}/confreec-${samp_name}-${ref_name}.txt"
  sed "s|strclf|${base_dir}/ref/chrseq-${ref_name}/len-chr.txt|" \
  < "${tplt_file}" | \
  sed "s|strcf|${base_dir}/ref/chrseq-${ref_name}|" | \
  sed "s|strmf|${base_dir}/map-sr/${samp_name}-${ref_name}-${bam_sfx}|" | \
  sed "s|strod|${samp_out_dir}|" | \
  sed "s|strmpf|${base_dir}/gem/mpt/${ref_name}.mappability|" \
  > "${my_config}"
    
  ### run
  freec -conf "${my_config}" &> "${samp_out_dir}/${samp_name}.log"
  ### format output file names
  file_gc=$(\ls "${samp_out_dir}/"*"GC_profile"*)
  win_size=$(basename "${file_gc}" | cut -d "." -f 2)
  mv -f "${file_gc}" \
  "${samp_out_dir}/${samp_name}-${ref_name}-gcp-${win_size}.txt"
  mv -f "${samp_out_dir}/"*"_ratio.txt" \
  "${samp_out_dir}/${samp_name}-ratio.txt"
  mv -f "${samp_out_dir}/"*"_info.txt" \
  "${samp_out_dir}/${samp_name}-info.txt"
  mv -f "${samp_out_dir}/"*"_CNVs" \
  "${samp_out_dir}/${samp_name}-cnvs.txt"
  mv -f "${samp_out_dir}/"*".cpn" \
  "${samp_out_dir}/${samp_name}-scn-profile.txt"
    
  ### calculate significance
  cat "${lib_dir}/sig-freec.R" | R --slave --args \
  "${samp_out_dir}/${samp_name}-cnvs.txt" \
  "${samp_out_dir}/${samp_name}-ratio.txt"
  rm "${samp_out_dir}/${samp_name}-cnvs.txt"
  mv -f "${samp_out_dir}/${samp_name}-cnvs.txt.p.value.txt" \
  "${samp_out_dir}/${samp_name}-cnvs.txt"
  
  ### clean plot output folder
  rm -f "${samp_out_dir}/"*"-ratio.txt"*"log2.pdf"
  
  ### plot
  all_chr=$(cut -f 2 "${base_dir}/ref/chrseq-${ref_name}/len-chr.txt" \
  | sed 's|chr||')
  n_chrom=$(echo "${all_chr}" | wc -w)
  cat "${lib_dir}/plot-freec.R" | R --slave --args "${n_ploidy}" \
  "${samp_out_dir}/${samp_name}-ratio.txt" "${n_chrom}" ${all_chr}
  
  ### rename the plots
  for ind_p in $(find "${samp_out_dir}/" -name "*ratio.txt*-log2.pdf"); do 
    new_name=$(echo "${ind_p}" | sed 's|-ratio.txt||')
    mv -f "${ind_p}" "${new_name}"
  done
  ) &
done

wait
