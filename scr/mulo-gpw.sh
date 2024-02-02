#!/bin/bash

## header ---------------------------------------------------------------------

### this script is the mulo-gpw runner

## user's settings ------------------------------------------------------------

ref_name="SK1"
n_ploidy="2"

### sample IDs where:
### - plus_samp is the (plus) selected sample in experimental design A or
###   the selected sample (plus or minus) in experimental design BH and BL
### - minu_samp is the (minus) selected sample in experimental design A or
###   the non-selected sample in experimental design BH and BL

### pm-1 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D2 D5 D8 D11"
### pm-2 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D5 D8 D11 D2"
### pm-3 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D8 D11 D2 D5"
### pm-4 (with exp_design="A")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D11 D2 D5 D8"
### pe-1 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D3 D6 D9 D12"
### pe-2 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D6 D9 D12 D3"
### pe-3 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D9 D12 D3 D6"
### pe-4 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D12 D3 D6 D9"
### me-1 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D3 D6 D9 D12"
### me-2 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D6 D9 D12 D3"
### me-3 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D9 D12 D3 D6"
### me-4 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D12 D3 D6 D9"
### pz-1 (with exp_design="BH")
# plus_samp="D1 D4 D7 D10"
# minu_samp="D13 D13 D13 D13"
### mz-1 (with exp_design="BL")
# plus_samp="D2 D5 D8 D11"
# minu_samp="D13 D13 D13 D13"
### ez-1 (with exp_design="BH", assuming "equal" samples have a plus phenotype)
plus_samp="D3 D6 D9 D12"
minu_samp="D13 D13 D13 D13"

### the experimental design ("exp_design" variable) can be "A"
### (for plus/minus phenotyped samples)
### or "BH" (for plus-evolved/non-evolved samples)
### or "BL" (for minus-evolved/non-evolved samples)
exp_design="BH"

## system's settings ----------------------------------------------------------

### check logs folder
if [[ ! -d "logs" ]]; then mkdir "logs"; fi

## clmnt ----------------------------------------------------------------------

### quality check
# bash fq-check.sh \
# > "logs/fq-check.out" 2> "logs/fq-check.err" &

### short-reads subshell
(
### reference indexing
bash index-ref.sh "${ref_name}" \
> "logs/index-ref.out" 2> "logs/index-ref.err"

### mapping
bash map-sr.sh "${ref_name}" \
> "logs/map-sr.out" 2> "logs/map-sr.err"
) &

### assemblies subshell
(
### align the parental genomes
bash parent-aligner.sh "${ref_name}" \
> "logs/parent-aligner.out" 2> "logs/parent-aligner.err"

### indexing (for vcf construction)
bash make-fai.sh \
> "logs/make-fai.out" 2> "logs/make-fai.err"

### transform nucmer data to vcf
Rscript nuc-vcf.R \
> "logs/nuc-vcf.out" 2> "logs/nuc-vcf.err"

### build the multisample vcf and the header for the fix multisample file
bash nuc-multivcf.sh \
> "logs/nuc-multivcf.out" 2> "logs/nuc-multivcf.err"

### fix the multisample vcf with nucmer alignment data
Rscript fix-gen.R \
> "logs/fix-gen.out" 2> "logs/fix-gen.err"
) &

wait

### extract private markers
Rscript private-markers.R
> "logs/private-markers.out" 2> "logs/private-markers.err"

## genotyping samples at marker positions
bash gen-markers.sh "${n_ploidy}" \
> "logs/gen-markers.out" 2> "logs/gen-markers.err"

### plot-stats.sh runs smooth-plot.R and BRM.R for each couple of samples 
### e.g. one plus and the corresponsding minus
bash plot-stats.sh "${plus_samp}" "${minu_samp}" "${exp_design}" \
> "logs/plot-stats.out" 2> "logs/plot-stats.err"

echo "Ah, pardon! Tarapio tapioco come se fosse antani, la supercazzola prematurata con dominus vobiscum blinda??"
echo "[Earl Conte Raffaello 'Lello' Mascetti]"
