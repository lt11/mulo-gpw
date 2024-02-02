## header ---------------------------------------------------------------------

### This script loads the multisample vcf file derived from nucmer alignments
### and, for each sample reported in the header, it extracts the loci
### that are private for one sample. Then it also extracts the loci where 
### the reference is the only genome bearing the REF allele.

rm(list = ls())
options(stringsAsFactors = F)
library(here)

## settings -------------------------------------------------------------------

### fixed settings
dirBase <- dirname(here())
### def off
# dirBase <- "/Users/Lorenzo/dev/dev-mulo-gpw"

### input and output folders
dirOut <- file.path(dirBase, "var-calls", "flt-markers")
unlink(dirOut, recursive = T)
dir.create(dirOut, recursive = T)
dirIn <- file.path(dirBase, "gen", "ms-par-vcf")

### output files prefix
outPref <- "private-"

### vcf data-frame fixed column names (no samples)
hdVcfFix <- c("Chrom_id", "Pos_bp", "Var_id", "Ref_allele", "Alt_allele",
              "Qual_val", "Filter_tag", "Info_tags", "Format_def")
### vcf header
headerVcf <- c("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")

## clmnt ----------------------------------------------------------------------

### read the vcf
fileMultis <- list.files(path = dirIn, pattern = "genfix.vcf$", full.names = T)
dfMultisVcf <- read.table(file = fileMultis)

## read the header of the vcf -------------------------------------------------

### build the command line to generate the header
fileHd <- paste0("hd-",
                 sub(pattern = "vcf",
                     replacement = "txt",
                     basename(fileMultis)))
pathFileHd <- file.path(dirname(fileMultis), fileHd)
makeHeader <- paste("grep ^#C",
                    fileMultis,
                    ">",
                    pathFileHd,
                    sep = " ")
system(makeHeader)
strVcf <- readLines(pathFileHd)
vtVcfHeader <- unlist(strsplit(strVcf, split = "\t"))

### build the command line to generate the header descriptor
fileHdDes <- c("hd-descriptor.txt")
pathFileHdDes <- file.path(dirname(fileMultis), fileHdDes)
makeHeaderDes <- paste("grep ^##",
                       fileMultis,
                       ">",
                       pathFileHdDes,
                       sep = " ")
system(makeHeaderDes)
strHdDescriptor <- readLines(pathFileHdDes)

### rename the columns
nFields <- length(vtVcfHeader)
nSamples <- nFields - 9
sampleFields <- vtVcfHeader[c(10:nFields)]
colnames(dfMultisVcf)[1:9] <- hdVcfFix
colnames(dfMultisVcf)[10:nFields] <- sampleFields

### keep only the GT subfield
strFormat <- unlist(strsplit(dfMultisVcf$Format_def[1], split = ":"))
indGT <- which(strFormat == "GT")
dfMultisVcf$Format_def <- rep("GT", nrow(dfMultisVcf))

for (indS in c(10:nFields)) {
  strGT <- sapply(strsplit(dfMultisVcf[, indS], split = ":"), "[[", indGT)
  ### indexes of missing genotypes
  indMissing <- which(strGT == ".")
  strGT[indMissing] <- "0"
  ### GT is numeric, this is very helpful to get the private
  ### since the correspond to, e.g. if we have 4 samples,
  ### to the sum of 4 or 0 
  ### (but missing GT must be converted from "." to 0)
  numGT <- as.numeric(strGT)
  dfMultisVcf[, indS] <- numGT
}
mtGT <- as.matrix(dfMultisVcf[, c(10:nFields)])
sumGT <- rowSums(mtGT)

## reference private loci -----------------------------------------------------

### indexes of the loci private to the reference
indRefPrivate <- which(sumGT == 4)
### make the vcf of the loci private to the reference
dfPrvRefVcf <- data.frame(dfMultisVcf[indRefPrivate, c(1:9)],
                          rep(0, length(indRefPrivate)))
### get the reference name
refFile <- basename(sub(pattern = "##reference=",
                        replacement = "",
                        x = grep(pattern = "##reference=",
                                 x = strHdDescriptor,
                                 value = T)))
refName <- sub("^([^-]*)-.*$", "\\1", refFile)
colnames(dfPrvRefVcf)[10] <- refName
### write the header descriptor (strHdDescriptor)
fileOutRef <- file.path(dirOut, paste0(outPref, refName, ".vcf"))
cat(strHdDescriptor, file = fileOutRef, sep = "\n")
### write the header (headerVcfRef)
headerVcfRef <- paste0(headerVcf, refName)
cat(headerVcfRef, file = fileOutRef, sep = "\n", append = T)
### write the single-sample vcf (without the data-frame column names)
write.table(x = dfPrvRefVcf, file = fileOutRef,
            append = T, quote = F, sep = "\t", col.names = F, row.names = F)

## samples private loci -------------------------------------------------------

### the weird initialisation of an empty list
lsPrivate <- vector(mode = "list", length = nSamples)
for (indG in c(1:nSamples)) {
  lsPrivate[[indG]] <- which(mtGT[, indG] == 1
                             & sumGT == 1)
  singleSampleVcf <- dfMultisVcf[lsPrivate[[indG]], c(1:9, c(9 + indG))]
  ### get who is the sample and make the name of the output file
  fileOutSampG <- file.path(dirOut, paste0(outPref, sampleFields[indG], ".vcf"))
  ### write the header descriptor (strHdDescriptor)
  cat(strHdDescriptor, file = fileOutSampG, sep = "\n")
  ### write the header (headerVcfSampG)
  headerVcfSampG <- paste0(headerVcf, sampleFields[indG])
  cat(headerVcfSampG, file = fileOutSampG, sep = "\n", append = T)
  ### write the single-sample vcf (without column names)
  write.table(x = singleSampleVcf, file = fileOutSampG,
              append = T, quote = F, sep = "\t", col.names = F, row.names = F)
}
