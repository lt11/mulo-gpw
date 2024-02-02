## header  --------------------------------------------------------------------

### the one that prepares the fasta files for GC-content calculation
 
## settings  ------------------------------------------------------------------

full_dir=$(cd $(dirname "${0}") && pwd)
base_dir=$(dirname "${full_dir}")
ref_name="${1}"

### working directory
ref_dir="${base_dir}/ref"

## clmnt  ---------------------------------------------------------------------

echo "Running the one that prepares the fasta files \
for GC-content calculation..."

for ind_ref in ${ref_name}; do
  out_dir="${base_dir}/ref/chrseq-${ind_ref}"
  if [[ -d "${out_dir}" ]]; then rm -rf "${out_dir}"; fi
  mkdir "${out_dir}"
  
  ref_fa=$(find "${ref_dir}" -name "${ind_ref}-genome.fa")
  ref_fai="${ref_fa}.fai"
  
  ### extract chromosomes IDs
  all_chr=$(grep "^>" "${ref_fa}" | sed 's|^.||' | cut -d" " -f 1)
  
  ### filter out small contigs and mitochondrial sequences
  my_chr=$(echo "${all_chr}" | \
  grep -v "mating_type_region" | \
  grep -v "chr_II_telomeric_gap" | \
  grep -v "chrMT" | \
  grep -v "MT" | \
  grep -v "mitochondrial")
    
  for ind_chr in ${my_chr}; do
    ### prepare *len-chr.txt
    grep -w "${ind_chr}" "${ref_fai}" | \
    awk 'BEGIN {FS="\t"} {OFS="\t"} {print NR, $1, $2}' \
    >> "${out_dir}/len-chr.txt"
    ### prepare the single-chromosome fasta files
    samtools faidx "${ref_fa}" "${ind_chr}" \
    > "${out_dir}/${ind_chr}.fa"
  done
done


