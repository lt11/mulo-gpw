## header ---------------------------------------------------------------------

### This script works on the merged files created by its runner 
### (e.g. D1-D2-SK1-mrg.txt.gz) which is stored in "/.../var-calls/mrg-markers". 
### For each pair (e.g. D1-SK1 and D2-SK1), it loads the corresponding vcf files 
### and filter them removing loci which are not genotyped on both files. 
### Then it removes also the loci with an ALT allele
### different from the allele reported in the vcf of the ancestor. 
### The latter is the vcf file which was filtered for e.g. SK1 private loci.
### A founder sample is then used to - sort of - normalise the AF: its AF 
### is subtracted from the AF of the samples (both the plus and the minus).
### The AF data are smoothed with a loess. 
### The smoothing and the AF data are finally plotted.
### It also prepares the data (in bsa format) for the statistical model of BRM.

rm(list = ls())
options(stringsAsFactors = F)
options(scipen = 999)
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
### dirMerged <- "/Users/Lorenzo/dev/dev-mulo-gpw/var-calls/mrg-markers/"
### fileMerged <- file.path(dirMerged, "D4-D5-DBVPG6044-mrg.txt.gz")

if (length(argsVal) == 2) {
  fileBedZoom <- argsVal[2]
  booZoom <- T
} else if (length(argsVal) > 2) {
  stop("Too many arguments provided!", call. = F)
}

### fixed settings
dirBase <- dirname(here())
### dev
### dirBase <- "/Users/Lorenzo/dev/dev-mulo-gpw"
spanVal <- 0.2

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
dirNucVcf <- file.path(dirBase, "gen", "ms-par-vcf")
### this folder is created in the plot-stats.sh
dirOutBrm <- file.path(dirBase, "allele-shift", "brm-input")
### this folder is created in the plot-stats.sh
dirDataAf <- file.path(dirBase, "allele-shift", "af")
### this folder is created in the plot-stats.sh
dirMergedMultiAl <- file.path(dirBase, "var-calls", "mrg-markers-multialleles")
### this folder is created in the plot-stats.sh
dirSingleSampleDifAlt <- file.path(dirBase, "var-calls", "gen-markers-diffalt")
### this folder is created in the plot-stats.sh
dirOutPlotDiff <- file.path(dirBase, "allele-shift", "plots", "af-diff")
### this folder is created in the plot-stats.sh
dirOutPlotFreq <- file.path(dirBase, "allele-shift", "plots", "af")

### in case we want to save the smoothed data 
### dirOutData <- file.path(dirBase, "allele-shift", "smoothed-data")
### dir.create(dirOutData, showWarnings = F, recursive = T)

### vcf file column names
hdVcfFix <- c("Chrom_id", "Pos_bp", "Var_id", "Ref_allele", "Alt_allele",
              "Qual_val", "Filter_tag", "Filter_tag", "Format_tag")

## clmnt ----------------------------------------------------------------------

### get the ID of the reference
dirRef <- file.path(dirBase, "ref")
refID <- sub(pattern = "^([^-]*)-.*$", replacement = "\\1",
             list.files(path = dirRef, pattern = "-genome.fa$"))

fileChromLen <- list.files(path = dirRef, pattern = "-genome.fa.fai$",
                           full.names = T)
dtChromLen <- fread(fileChromLen, sep = "\t")[, 1:2]
colnames(dtChromLen) <- c("Chrom_id", "Len_bp")

pairID <- sub("^([^-]*-[^-]*-[^-]*)-.*$", "\\1", basename(fileMerged))
cat("Processing sample pair: ", pairID, "\n", sep = "")
### get the ID of ancestor (we are working on its private alleles)
ancPrvAlle <- sub("^([^-]*-[^-]*)-(.*$)", "\\2", pairID)

### load the merged file (e.g. D1-D2-SK1-mrg.txt.gz)
dtMrgVcf <- fread(file = fileMerged, sep = "\t", header = T)
### fix the header
sampleIDs <- basename(colnames(dtMrgVcf)[c(10, 11)])
sampleIDs <- sub("^([^-]*)-.*$", "\\1", sampleIDs)
hdVcfFull <- c(hdVcfFix, sampleIDs)
colnames(dtMrgVcf) <- hdVcfFull

## filter #1: keep loci genotyped in both samples -----------------------------

cat("Filtering out loci which were not genotyped in both samples...",
    "\n", sep = "")
nAncLoci <- nrow(dtMrgVcf)
cat("Number of loci in the merged file: ", nAncLoci, "\n", sep = "")
nCol <- ncol(dtMrgVcf)
nSamples <- nCol - 9
### number of samples with the variant (locus by locus)
nSampVar <- apply(dtMrgVcf[, 10:nCol], 1, PresVar)
diffSamp <- nSamples - nSampVar
indBonGen <- which(diffSamp == 0)
indMissingGen <- which(diffSamp != 0)
### fraction of variants retained
fracRet <- length(indBonGen) / nAncLoci
cat("Fraction of loci retained: ", fracRet, "\n", sep = "")
### and finally filtering out
if (length(indMissingGen) != 0) {
  dtMrgVcf <- dtMrgVcf[-indMissingGen, ]
}

## filter #2: check if there is any multiallelic locus in dtMrgVcf ------------

indMulti <- which(nchar(dtMrgVcf[["Alt_allele"]]) != 1)
if (length(indMulti) != 0) {
  dtMrgVcfMultiAllele <- dtMrgVcf[indMulti, ]
  ### write to disk
  fileMultiAl <- paste0(pairID, "-mrg-multialleles.txt")
  pathMultiAl <- file.path(dirMergedMultiAl, fileMultiAl)
  fwrite(dtMrgVcfMultiAllele, file = pathMultiAl, quote = F, sep = "\t",
         row.names = F, col.names = T)
  ### filter out the multiallelic loci
  dtMrgVcf <- dtMrgVcf[-indMulti, ]
}

## filter #3: keep only the loci reported in dtMrgVcf -------------------------

### load the single-sample vcf files
plusID <- sampleIDs[1]
minuID <- sampleIDs[2]
plusPatt <- paste0(plusID, "-", ancPrvAlle, "\\.vcf\\.gz$")
minuPatt <- paste0(minuID, "-", ancPrvAlle, "\\.vcf\\.gz$")
filePlus <- list.files(path = dirSingleSample,
                       pattern = plusPatt,
                       full.names = T)
fileMinu <- list.files(path = dirSingleSample,
                       pattern = minuPatt,
                       full.names = T)
dtVcfPlus <- fread(file = filePlus, sep = "\t", header = T)
### for unknown reasons fread refuses to read the vcf file
### created by awk rather than a standard tool
### (although it is perfectly formatted);
### here this can happen when we use the founder data 
### instead of minus samples (i.e. when we run in evolved/founder mode
### instead of the plus/minus mode)
rileccaMelo <- read.table(file = fileMinu, header = F, sep = "\t")
dtVcfMinu <- data.table(rileccaMelo)
rm(rileccaMelo)
### fix the headers
hdVcfPlus <- c(hdVcfFix, plusID)
colnames(dtVcfPlus) <- hdVcfPlus
hdVcfMinu <- c(hdVcfFix, minuID)
colnames(dtVcfMinu) <- hdVcfMinu

### load the private loci (dtAncPrv) of an ancestor
prvLociPatt <- paste0("private-", ancPrvAlle, "\\.vcf$")
fileAncFlt <- list.files(path = dirAncFlt, pattern = prvLociPatt,
                         full.names = T)
### for unknown reasons fread refuses to read the vcf file
### created by awk rather than a standard tool
### (although it is perfectly formatted)
leccaMelo <- read.table(file = fileAncFlt, header = F, sep = "\t")
dtAncPrv <- data.table(leccaMelo)
### fix the header
hdVcfAnc <- c(hdVcfFix, ancPrvAlle)
colnames(dtAncPrv) <- hdVcfAnc

### filtering
# wanna test it?
# aa <- c("a","b","c","d","e","f")
# bb <- c("b","d","f")
# match(setdiff(aa, bb), aa)
strMrgVcf <- paste0(dtMrgVcf[["Chrom_id"]], dtMrgVcf[["Pos_bp"]])
strVcfPlus <- paste0(dtVcfPlus[["Chrom_id"]], dtVcfPlus[["Pos_bp"]])
strVcfMinu <- paste0(dtVcfMinu[["Chrom_id"]], dtVcfMinu[["Pos_bp"]])
strAncPrv <- paste0(dtAncPrv[["Chrom_id"]], dtAncPrv[["Pos_bp"]])
indBadVcfPlus <- match(setdiff(strVcfPlus, strMrgVcf), strVcfPlus)
indBadVcfMinu <- match(setdiff(strVcfMinu, strMrgVcf), strVcfMinu)
indBadAncPrv <- match(setdiff(strAncPrv, strMrgVcf), strAncPrv)
### we do not keep track of indBad* since we have already saved the
### multiallelic loci and here we have considered those which are not
### genotyped in all samples due to low coverage
if (length(indBadVcfPlus) != 0) {
  dtVcfPlus <- dtVcfPlus[-indBadVcfPlus, ]
}
if (length(indBadVcfMinu) != 0) {
  dtVcfMinu <- dtVcfMinu[-indBadVcfMinu, ]
}
if (length(indBadAncPrv) != 0) {
  dtAncPrv <- dtAncPrv[-indBadAncPrv, ]
}

## filter #4: check if ALTs in single-sample vcf files match with dtAncPrv ----

### since we use a bit of non-base R better have a quick check: 
### are the vcf files well sorted?
checkPlusFile <- which(dtAncPrv[["Pos_bp"]] != dtVcfPlus[["Pos_bp"]])
checkMinuFile <- which(dtAncPrv[["Pos_bp"]] != dtVcfMinu[["Pos_bp"]])
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

### take care: we need to check only the ALT alleles that are not "."; e.g. if
### the minus sample has only REF alleles the ALT column will report "."
### which is totally fine (it means that the ALT allele present in the ancestor,
### and enriched or fixed in the plus population, is partially or completely 
### missing in the minus population)

### check if the ALT allele in dtVcfPlus matches with dtAncPrv
### and in case keep those which do not match
indBadAltPlus <- which(dtAncPrv[["Alt_allele"]] != dtVcfPlus[["Alt_allele"]]
                       & dtVcfPlus[["Alt_allele"]] != ".")
if (length(indBadAltPlus) != 0) {
  dtVcfPlusDiffAlt <- dtVcfPlus[indBadAltPlus, ]
  fileDifAlt <- paste0(plusID, "-", ancPrvAlle, "-diffalt.txt")
  pathOut <- file.path(dirSingleSampleDifAlt, fileDifAlt)
  fwrite(dtVcfPlusDiffAlt, file = pathOut, quote = F, sep = "\t",
         row.names = F, col.names = T)
}

### check if the ALT allele in dtVcfMinu matches with dtAncPrv
### and in case put apart those which do not match
indBadAltMinu <- which(dtAncPrv[["Alt_allele"]] != dtVcfMinu[["Alt_allele"]]
                       & dtVcfMinu[["Alt_allele"]] != ".")
if (length(indBadAltMinu) != 0) {
  dtVcfMinuDiffAlt <- dtVcfMinu[indBadAltMinu, ]
  fileDifAlt <- paste0(minuID, "-", ancPrvAlle, "-diffalt.txt")
  pathOut <- file.path(dirSingleSampleDifAlt, fileDifAlt)
  fwrite(dtVcfMinuDiffAlt, file = pathOut, quote = F, sep = "\t",
         row.names = F, col.names = T)
}

### set operations
### union elements are unique (no duplicates)
indBadAlt <- union(indBadAltPlus, indBadAltMinu)
if (length(indBadAlt) != 0) {
  dtVcfPlus <- dtVcfPlus[-indBadAlt, ]
  dtVcfMinu <- dtVcfMinu[-indBadAlt, ]
  dtAncPrv <- dtAncPrv[-indBadAlt, ]
}

## TODO: normalise the AF with the AF values of the founder -------------------

### è il caso di normalizzare in qualche modo coi dati delle AF dell'ancestor
### (che non abbiamo genotipizzato)?
### sì, ma scrivilo in modo che possa essere cancellato con un colpo di CANC
### va cambiato anche la riga del wrapper 
### bash plot-stats.sh "${plus_samp}" "${minu_samp}" "${exp_design}"
### perché gli servono gli ID dei founder

### in caso rammentati che vanno usati solo i loci sopravvissuti ai filtracci
### in uno a scelta tra dtVcfPlus, dtVcfMinu, o dtAncPrv (hanno tutti gli
### stessi loci) in comune con quelli presenti nel vcf del founder

### poi c'è da controllare che alleli ALT (quelli diversi da ".") nel founder
### siano uguali agli ALT di dtAncPrv (chiamiamoli indBadFound)

### e infine vanno tolti eventuali alleli indBadFound da, ovviamente il 
### founder, ma anche da dtVcfPlus, dtVcfMinu, e dtAncPrv

## calculate AF from single-sample vcf, using the DP4 values ------------------

### we use the single-sample files since the DP4 tag is in the INFO field
### and we trust it more than the DP tag (which in principle would be easier 
### to handle since it is reported in the FORMAT field)
countPlus <- ExtractDP4(dtVcfPlus)
countMinu <- ExtractDP4(dtVcfMinu)

## data for BRM ---------------------------------------------------------------

### of course, dtVcfPlus$Chrom_id == dtVcfMinu$Chrom_id
dtWholeBsa <- data.table(dtVcfPlus$Chrom_id, dtVcfPlus$Pos_bp,
                         countPlus$countRefDP4, countPlus$countAltDP4,
                         countMinu$countRefDP4, countMinu$countAltDP4)
colnames(dtWholeBsa) <- c("Chrom_id", "Pos_bp",
                          "Plus_ref_count", "Plus_alt_count",
                          "Minus_ref_count", "Minus_alt_count")
dtBsa <- dtWholeBsa[Chrom_id != "chrMT"]
### save the bsa data
fileBsaOut <- paste0(pairID, "-bsa.txt")
pathBsaOut <- file.path(dirOutBrm, fileBsaOut)
fwrite(x = dtBsa, file = pathBsaOut, sep = "\t", quote = F,
       row.names = F, col.names = F)

## data for ggplot ------------------------------------------------------------

refFracPlus <- countPlus$countRefDP4 / 
  (countPlus$countRefDP4 + countPlus$countAltDP4)
refFracMinu <- countMinu$countRefDP4 / 
  (countMinu$countRefDP4 + countMinu$countAltDP4)

dtPlot <- data.table(dtVcfPlus$Chrom_id,
                     dtVcfPlus$Pos_bp,
                     refFracPlus,
                     refFracMinu)
colnames(dtPlot)[1:4] <- c("Chr_id", "Pos_kbp",
                           "Ref_frac_plus", "Ref_frac_minus")

### switch to kbp
dtPlot$Pos_kbp <- dtPlot$Pos_kbp / 1000
fileDataOut <- file.path(dirDataAf, paste0(pairID, "-af.RData"))
save(dtPlot, file = fileDataOut)

## plotting -------------------------------------------------------------------

for (indC in unique(dtPlot$Chr_id)) {
  indL <- which(dtChromLen$Chrom_id == indC)
  coordStart <- 0.001
  coordEnd <- dtChromLen$Len_bp[indL] / 1000
  ### AF plot
  fileFreqName <- paste0(plusID, "-", minuID, "-",
                         ancPrvAlle, "-", indC, ".pdf")
  fileFreqOut <- file.path(dirOutPlotFreq, fileFreqName)
  dtPlotChr <- dtPlot[indC, on = "Chr_id"]
  plotTolo <- ggplot(dtPlotChr) +
    theme_minimal() +
    theme(plot.margin = unit(c(1, 2, 1, 1), "cm"),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 22)) +
    coord_cartesian(ylim = c(.0, 1.), xlim = c(coordStart, coordEnd)) +
    labs(title = "", subtitle = "",
         x = paste0("Position (", refID, " ", indC,") [kbp]"), y = "AF",
         size = 28) +
    ### plus sample
    geom_point(aes(Pos_kbp, dtPlotChr[[3]]), colour = "red",
               size = 0.5, alpha = 0.3) +
    ### stat_smooth may generate warnings if there are a few loci
    stat_smooth(aes(Pos_kbp, dtPlotChr[[3]]), method = "loess",
                span = spanVal, colour = "red", alpha = 0.6,
                formula = y ~ x,
                n = 80) +
    ### minus sample
    geom_point(aes(Pos_kbp, dtPlotChr[[4]]), colour = "blue",
               size = 0.5, alpha = 0.3) +
    stat_smooth(aes(Pos_kbp, dtPlotChr[[4]]), method = "loess",
                span = spanVal, colour = "blue", alpha = 0.6,
                formula = y ~ x,
                n = 80)
  pdf(file = fileFreqOut)
  print(plotTolo)
  dev.off()
  
  ### AF difference plot
  fileDiffName <- paste0("afd-", plusID, "-", minuID, "-",
                         ancPrvAlle, "-", indC, ".pdf")
  fileDiffOut <- file.path(dirOutPlotDiff, fileDiffName)
  plotTolo <- ggplot(dtPlotChr) +
    theme_minimal() +
    theme(plot.margin = unit(c(1, 2, 1, 1), "cm"),
          axis.text = element_text(size = 18),
          axis.title = element_text(size = 22)) +
    coord_cartesian(ylim = c(-0.5, +0.5), xlim = c(coordStart, coordEnd)) +
    labs(title = "", subtitle = "",
         x = paste0("Position (", refID, " ", indC,") [kbp]"),
         y = "AFD (plus - minus)",
         size = 28) +
    ### difference
    geom_point(aes(Pos_kbp, dtPlotChr[[3]] - dtPlotChr[[4]]), colour = "black",
               size = 0.5, alpha = 0.3) +
    stat_smooth(aes(Pos_kbp, dtPlotChr[[3]] - dtPlotChr[[4]]), method = "loess",
                span = spanVal, colour = "black", alpha = 0.6,
                formula = y ~ x,
                n = 80)
  pdf(file = fileDiffOut)
  print(plotTolo)
  dev.off()
}
