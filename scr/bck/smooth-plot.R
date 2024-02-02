## header ---------------------------------------------------------------------

### This script works on the merged files created by its runner 
### (e.g. D1-D2-mrg.txt.gz) which is stored in "/.../var-calls/mrg-markers".
### For each pair (e.g. D1 and D2),
### it loads the corresponding vcf files and filter 
### them removing loci which are not
### genotyped on both files. Then it removes also the loci with an ALT allele
### different from the allele reported in the vcf of the ancestor. The latter
### is the vcf file which was quality-filtered. Then the vcf files are used
### to calculate the allele frequency (AF) at each locus. The AF data are
### smoothed with a loess. The smoothing and the AF data are finally plotted.
### It also prepares the data (in bsa format) for the statistical model of BRM.

rm(list = ls())
options(stringsAsFactors = F)
library(stringr)
library(ggplot2)
library(here)
library(data.table)

## function(s) ----------------------------------------------------------------

### A note about PresVar:
### if a variant is present in a vcf file it is reported 
### as e.g. "1:114,13:4:0:0,0:1,3:1,3"
### while if it not present the locus reports 
### ".:.:.:.:.:.:."
### so if we want to know if a locus has been genotyped
### we can just look for digits (i.e. 0, 1, 2, …, 9) 
### in the "SAMPLE" string of the vcf

PresVar <- function(x) length(grep("[[:digit:]]+", x))
### It searches for digits in a vector of strings.
### 
### Arguments:
### (1) a vector of string
###
### Returns:
### (1) an integer number (the number of elements with at least one digits)

ExtractDP4 <- function(x) {
### It reads a vcf (with the DP4 tag embedded in the INFO field)
### in a data table and returns a data table with two columns
### with the DP4 of the REF and the ALT alleles.
### One more example: if a INFO entry contains "DP4=92,81,35,44"
### then strTagDP4 is "67,81,45,47".
### 
### arguments:
### (1) a vcf in a data table
###
### returns:
### (1) a data table with two columns with the DP4 
### of the REF and the ALT alleles
  strTagDP4 <- sub(pattern = "^.*DP4=([^;]*);.*$",
                   replacement = "\\1",
                   x = x[[8]]) ### the 8th column of the data table
  refForward <- as.numeric(sapply(strsplit(strTagDP4, split = ","), "[[", 1))
  refReverse <- as.numeric(sapply(strsplit(strTagDP4, split = ","), "[[", 2))
  altForward <- as.numeric(sapply(strsplit(strTagDP4, split = ","), "[[", 3))
  altReverse <- as.numeric(sapply(strsplit(strTagDP4, split = ","), "[[", 4))
  countRefDP4 <- refForward + refReverse
  countAltDP4 <- altForward + altReverse
  dtDP4 <- data.table(countRefDP4, countAltDP4)
  return(dtDP4)
}

## settings -------------------------------------------------------------------

### default values
booZoom <- F
fileBedZoom <- "dummy"

### arguments
argsVal <- commandArgs(trailingOnly = T)
fileMerged <- argsVal[1]
### dev
### fileMerged <- "/Users/Lorenzo/Desktop/dev-mulo-ea/var-calls/mrg-markers/D3-D13-mrg.txt.gz"

if (length(argsVal) == 2) {
  fileBedZoom <- argsVal[2]
  booZoom <- T
} else if (length(argsVal) > 2) {
  stop("Too many arguments provided, only one argument is allowed!", call. = F)
}

### fixed settings
dirBase <- dirname(here())
### dev
### dirBase <- "/Users/Lorenzo/Desktop/dev-mulo-ea"
spanVal <- 0.3

### optional zoom-in in regions defined in an (optional) input bed file
zoomData <- c()
if (file.exists(fileBedZoom)) {
  zoomData <- read.table(file = fileBedZoom, header = F, sep = "\t")
  zoomDataBp <- zoomData
  zoomDataBp[, 2] <- as.character(zoomDataBp[, 2])
  zoomDataBp[, 3] <- as.character(zoomDataBp[, 3])
  zoomData[, 2] <- zoomData[, 2] / 1000
  zoomData[, 3] <- zoomData[, 3] / 1000
}

### input and output folders
dirMerged <- file.path(dirBase, "var-calls", "mrg-markers")
dirSingleSample <- file.path(dirBase, "var-calls", "gen-markers")
dirAncFlt <- file.path(dirBase, "var-calls", "flt-markers")
dirOutPlotFreq <- file.path(dirBase, "allele-shift", "plots", "af")
dirOutPlotDiff <- file.path(dirBase, "allele-shift", "plots", "af-diff")
dirOutBrm <- file.path(dirBase, "allele-shift", "brm-input")
dir.create(dirOutPlotDiff, showWarnings = F, recursive = T)
dir.create(dirOutPlotFreq, showWarnings = F, recursive = T)
dir.create(dirOutBrm, showWarnings = F, recursive = T)
### in case we want to save the smoothed data 
### dirOutData <- file.path(dirBase, "allele-shift", "data")
### dir.create(dirOutData, showWarnings = F, recursive = T)

### vcf file column names
hdVcfFix <- c("chrom_id", "pos_bp", "var_id", "ref_allele", "alt_allele",
              "qual_val", "filter_tag", "info_tags", "format_def")

## clmnt ----------------------------------------------------------------------

### get the ID of the reference
dirRef <- file.path(dirBase, "ref")
refID <- sub(pattern = "^([^-]*)-.*$", replacement = "\\1",
             list.files(path = dirRef, pattern = "-genome.fa$"))

pairID <- sub("^([^-]*-[^-]*)-.*$", "\\1", basename(fileMerged))
cat("Processing sample pair: ", pairID, "\n", sep = "")
### get the chromosome IDs from the merged files (e.g. D1-D2-mrg.txt.gz)
mergedVcf <- fread(file = fileMerged, sep = "\t", header = T)
### fix the header
sampleIDs <- colnames(mergedVcf)[c(10, 11)]
sampleIDs <- sub("^([^-]*)-.*$", "\\1", sampleIDs)
hdVcfFull <- c(hdVcfFix, sampleIDs)
colnames(mergedVcf) <- hdVcfFull
  
## filters (for markers) ----------------------------------------------------

### filter #1: we thrash the loci which were not genotyped in both samples
cat("Filtering out loci which were not genotyped in all samples...",
    "\n", sep = "")
nAncLoci <- nrow(mergedVcf)
cat("Number of ancestral loci: ", nAncLoci, "\n", sep = "")
nCol <- ncol(mergedVcf)
nSamples <- nCol - 9
### number of samples with the variant (locus by locus)
nSampVar <- apply(mergedVcf[, 10:nCol], 1, PresVar)
diffSamp <- nSamples - nSampVar
indBonGen <- which(diffSamp == 0)
indMissingGen <- which(diffSamp != 0)
### fraction of variants retained
fracRet <- length(indBonGen) / nAncLoci
cat("Fraction of loci retained: ", fracRet, "\n", sep = "")
  
### filter #2
### than we control the alternative allele in the single-sample vcf files
### (vs e.g. D13-SK1-flt.vcf.gz)
### since here we still have the same number of variants in all the files
cat("Filtering out loci with ALT allele different from the ancestor...",
    "\n", sep = "")
  
### for unknown reasons fread refuses to read the vcf file
### created by awk rather than a standard tool
### (although it is perfectly formatted)
fileAncFlt <- list.files(path = dirAncFlt, pattern = "\\.vcf\\.gz$",
                         full.names = T)
leccaMelo <- read.table(file = fileAncFlt, header = F, sep = "\t")
vcfAncFlt <- data.table(leccaMelo)
### fix the header
ancID <- sub("^([^-]*)-.*$", "\\1", basename(fileAncFlt))
hdVcfAnc <- c(hdVcfFix, ancID)
colnames(vcfAncFlt) <- hdVcfAnc
  
### here we load the single-sample files
plusID <- sampleIDs[1]
minuID <- sampleIDs[2]
filePlus <- list.files(path = dirSingleSample,
                       pattern = paste0(plusID, "\\.vcf\\.gz$"),
                       full.names = T)
fileMinu <- list.files(path = dirSingleSample,
                       pattern = paste0(minuID, "\\.vcf\\.gz$"),
                       full.names = T)
vcfPlus <- fread(file = filePlus, sep = "\t", header = T)
### for unknown reasons fread refuses to read the vcf file
### created by awk rather than a standard tool
### (although it is perfectly formatted);
### here this can happen when we use the ancestral data 
### instead of minus samples (i.e. when we run in evolved/ancestral mode
### instead of the plus/minus mode)
rileccaMelo <- read.table(file = fileMinu, header = F, sep = "\t")
vcfMinu <- data.table(rileccaMelo)
### fix the headers
hdVcfPlus <- c(hdVcfFix, plusID)
colnames(vcfPlus) <- hdVcfPlus
hdVcfMinu <- c(hdVcfFix, minuID)
colnames(vcfMinu) <- hdVcfMinu
  
### a quick check (are the vcf files well sorted?)
checkPlusFile <- which(vcfAncFlt$pos_bp != vcfPlus$pos_bp)
checkMinuFile <- which(vcfAncFlt$pos_bp != vcfMinu$pos_bp)
if (length(checkPlusFile) != 0) {
  stop("It seems that the vcf file of the 'plus' sample 
       is not sorted as the ancestral one",
       call. = F)
}
if (length(checkMinuFile) != 0) {
  stop("It seems that the vcf file of the 'minus' sample 
       is not sorted as the ancestral one",
       call. = F)
}
  
indBadAltPlus <- which(vcfAncFlt$alt_allele != vcfPlus$alt_allele)
nBadAltPlus <- length(indBadAltPlus)
indBadAltMinu <- which(vcfAncFlt$alt_allele != vcfMinu$alt_allele)
nBadAltMinu <- length(indBadAltMinu)
  
### set operations
### union elements are unique (no duplicates)
unionBadAltLoci <- union(indBadAltMinu, indBadAltPlus)
nUnionBadAltLoci <- length(unionBadAltLoci)
onlyMinuBadAltLoci <- setdiff(indBadAltMinu, indBadAltPlus)
nOnlyMinuBadAltLoci <- length(onlyMinuBadAltLoci)
onlyPlusBadAltLoci <- setdiff(indBadAltPlus, indBadAltMinu)
nOnlyPlusBadAltLoci <- length(onlyPlusBadAltLoci)
sharedBadAltLoci <- intersect(indBadAltMinu, indBadAltPlus)
nSharedBadAltLoci <- length(sharedBadAltLoci)
  
cat("Number [%] of total 'bad-ALT' loci (both plus and minus samples): ",
    nUnionBadAltLoci,
    " [",
    format(round(nUnionBadAltLoci/nAncLoci*100, 1), nsmall = 1),
    "]",
    "\n", sep = "")
cat("Number [%] of minus-only 'bad-ALT' loci: ",
    nOnlyMinuBadAltLoci,
    " [",
    format(round(nOnlyMinuBadAltLoci/nAncLoci*100, 1), nsmall = 1),
    "]", "\n", sep = "")
cat("Number [%] of plus-only 'bad-ALT' loci: ",
    nOnlyPlusBadAltLoci,
    " [",
    format(round(nOnlyPlusBadAltLoci/nAncLoci*100, 1), nsmall = 1),
    "]",
    "\n", sep = "")
cat("Number [%] of shared (among plus and minus samples) 'bad-ALT' loci: ",
    nSharedBadAltLoci,
    " [",
    format(round(nSharedBadAltLoci/nAncLoci*100, 1), nsmall = 1),
    "]",
    "\n", sep = "")
### finally we remove the bad loci
indThrash <- unique(c(indMissingGen, unionBadAltLoci))
  
if (length(indThrash) != 0) {
  vcfMinu <- vcfMinu[-indThrash, ]
  vcfPlus <- vcfPlus[-indThrash, ]
  }

## calculation of AF from single-sample vcf files ---------------------------

### we use the single-sample files since the DP4 tag is the INFO field
### and we trust it more than the DP tag (which in principle would be easier 
### to handle since it is reported in the FORMAT field)

countPlus <- ExtractDP4(vcfPlus)
countMinu <- ExtractDP4(vcfMinu)

## data for BRM -------------------------------------------------------------

### of course, vcfPlus$chrom_id == vcfMinu$chrom_id
dtWholeBsa <- data.table(vcfPlus$chrom_id, vcfPlus$pos_bp,
                         countPlus$countRefDP4, countPlus$countAltDP4,
                         countMinu$countRefDP4, countMinu$countAltDP4)
colnames(dtWholeBsa) <- c("chrom_id", "pos_bp",
                          "plus_ref_count", "plus_alt_count",
                          "minus_ref_count", "minus_alt_count")
dtBsa <- dtWholeBsa[chrom_id != "chrMT"]
### save the bsa data
fileBsaOut <- paste0(pairID, "-bsa.txt")
pathBasOut <- file.path(dirOutBrm, fileBsaOut)
fwrite(x = dtBsa, file = pathBasOut, sep = "\t", quote = F,
       row.names = F, col.names = F)

## data for ggplot ----------------------------------------------------------

refFracPlus <- countPlus$countRefDP4 / 
  (countPlus$countRefDP4 + countPlus$countAltDP4)
refFracMinu <- countMinu$countAltDP4 / 
  (countMinu$countRefDP4 + countMinu$countAltDP4)

dtPlot <- data.table(vcfPlus$chrom_id,
                     vcfPlus$pos_bp,
                     refFracPlus,
                     refFracMinu)
colnames(dtPlot)[1:2] <- c("chr_id", "pos_kb")

### switch to kb
dtPlot$pos_kb <- dtPlot$pos_kb / 1000

## plotting -----------------------------------------------------------------

for (indC in unique(dtPlot$chr_id)) {
  ### AF plot
  fileFreqName <- paste0(plusID, "-", minuID, "-", indC, ".pdf")
  fileFreqOut <- file.path(dirOutPlotFreq, fileFreqName)
  dtPlotChr <- dtPlot[indC, on = "chr_id"]
  plotTolo <- ggplot(dtPlotChr) +
    theme_minimal() +
    theme(plot.margin = unit(c(1, 2, 1, 1), "cm"),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 22)) +
    coord_cartesian(ylim = c(.0, 1.)) +
    labs(title = "", subtitle = "",
         x = paste0("Position (",refID, " ", indC,") [kb]"), y = "AF",
         size = 28) +
    ### plus sample
    geom_point(aes(pos_kb, dtPlotChr[[3]]), colour = "red",
               size = 0.5, alpha = 0.3) +
    stat_smooth(aes(pos_kb, dtPlotChr[[3]]), method = "loess",
                span = spanVal, colour = "red", alpha = 0.6,
                formula = y ~ x,
                n = 80) +
    ### minus sample
    geom_point(aes(pos_kb, dtPlotChr[[4]]), colour = "black",
               size = 0.5, alpha = 0.3) +
    stat_smooth(aes(pos_kb, dtPlotChr[[4]]), method = "loess",
                span = spanVal, colour = "black", alpha = 0.6,
                formula = y ~ x,
                n = 80)
  pdf(file = fileFreqOut)
  print(plotTolo)
  dev.off()
  
  ### AF difference plot
  fileDiffName <- paste0("afd-", plusID, "-", minuID, "-", indC, ".pdf")
  fileDiffOut <- file.path(dirOutPlotDiff, fileDiffName)
  plotTolo <- ggplot(dtPlotChr) +
    theme_minimal() +
    theme(plot.margin = unit(c(1, 2, 1, 1), "cm"),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 22)) +
    coord_cartesian(ylim = c(-1., +1.)) +
    labs(title = "", subtitle = "",
         x = paste0("Position (",refID, " ", indC,") [kb]"),
         y = "AFD (plus - minus)",
         size = 28) +
    ### difference
    geom_point(aes(pos_kb, dtPlotChr[[3]] - dtPlotChr[[4]]), colour = "black",
               size = 0.5, alpha = 0.3) +
    stat_smooth(aes(pos_kb, dtPlotChr[[3]] - dtPlotChr[[4]]), method = "loess",
                span = spanVal, colour = "black", alpha = 0.6,
                formula = y ~ x,
                n = 80)
  pdf(file = fileDiffOut)
  print(plotTolo)
  dev.off()
  
}


