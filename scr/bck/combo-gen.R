## header ---------------------------------------------------------------------

options(stringsAsFactors = F)
library(here)

## settings -------------------------------------------------------------------

### arguments
argsVal <- commandArgs(trailingOnly = T)
plusSamples <- argsVal[1]
minuSamples <- argsVal[2]

### editable settings
plusThreshold <- 0.95
minuThreshold <- 0.05

### fixed settings
baseDir <-  dirname(here())
cat("Checking for positively selected markers in samples: ",
    plusSamples, "\n", sep = "")
cat("...AND...", "\n")
cat("Checking for negatively selected markers in samples: ",
    minuSamples, "\n", sep = "")

vcfDir <- file.path(baseDir, "var-calls" ,"gen-markers")
vcfHeader <- c("chr_str", "pos_num", "id_str", "ref_str", "alt_str",
               "qual_num", "filter_str", "info_str", "form_str", "samp_str")
outDir <- file.path(baseDir, "var-calls", "sel-markers")
dir.create(path = outDir, showWarnings = F)
outFile <- file.path(outDir, "candidates.txt")

## clmnt ----------------------------------------------------------------------

### split sample
plusSplit <- unlist(strsplit(x = plusSamples, split = " "))
minuSplit <- unlist(strsplit(x = minuSamples, split = " "))

### load all vcf files
allVcf <- list.files(path = vcfDir, full.names = T, pattern = "\\.vcf\\.gz$")

### read on vcf to get the number of markers
oneVcf <- read.table(file = allVcf[1], header = F)
nMarkers <- nrow(oneVcf)

### a matrix with plus samples BAF values, and one for the minus samples
bafPlusMx <- matrix(data = NA, nrow = nMarkers, ncol = length(plusSplit))
bafMinuMx <- matrix(data = NA, nrow = nMarkers, ncol = length(minuSplit))
### and the corresponding boolean matrices
bafBooPlusMx <- matrix(data = 0, nrow = nMarkers, ncol = length(plusSplit))
bafBooMinuMx <- matrix(data = 0, nrow = nMarkers, ncol = length(minuSplit))

counterP <- 1
for (indS in plusSplit) {
  myFile <- grep(pattern = indS, x = allVcf, value = T)
  myVcf <- read.table(file = myFile, header = F)
  colnames(myVcf) <- vcfHeader
  dp4Tag <- sub(pattern = "^.*(DP4=[^;]*).*$",
                replacement = "\\1",
                x = myVcf$info_str)
  dp4Val <- sub(pattern = "DP4=", replacement = "", dp4Tag)
  nForRef <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\1", x = dp4Val))
  nRevRef <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\2", x = dp4Val))
  nForAlt <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\3", x = dp4Val))
  nRevAlt <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\4", x = dp4Val))
  bafSamp <- (nForAlt + nRevAlt) / (nForRef + nRevRef + nForAlt + nRevAlt)
  bafPlusMx[, counterP] <- bafSamp
  ### make the boolean
  indBon <- which(bafSamp > plusThreshold)
  bafBooPlusMx[indBon, counterP] <- 1
  counterP <- counterP + 1
}

counterM <- 1
for (indS in minuSplit) {
  myFile <- grep(pattern = indS, x = allVcf, value = T)
  myVcf <- read.table(file = myFile, header = F)
  colnames(myVcf) <- vcfHeader
  dp4Tag <- sub(pattern = "^.*(DP4=[^;]*).*$",
                replacement = "\\1",
                x = myVcf$info_str)
  dp4Val <- sub(pattern = "DP4=", replacement = "", dp4Tag)
  nForRef <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\1", x = dp4Val))
  nRevRef <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\2", x = dp4Val))
  nForAlt <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\3", x = dp4Val))
  nRevAlt <- as.numeric(sub(pattern = "^([^,]*),([^,]*),([^,]*),([^,]*)$",
                            replacement = "\\4", x = dp4Val))
  bafSamp <- (nForAlt + nRevAlt) / (nForRef + nRevRef + nForAlt + nRevAlt)
  bafMinuMx[, counterM] <- bafSamp
  ### make the boolean
  indBon <- which(bafSamp < minuThreshold)
  bafBooMinuMx[indBon, counterM] <- 1
  counterM <- counterM + 1
}

totalMx <- cbind(bafBooPlusMx, bafBooMinuMx)
sumMx <- rowSums(totalMx)
indBingo <- which(sumMx == length(c(plusSplit, minuSplit)))

onePlusVcf <- grep(pattern = plusSplit[1], x = allVcf, value = T)
onePlusData <- read.table(file = onePlusVcf, header = F)
colnames(onePlusData) <- vcfHeader
outCoor <- data.frame(onePlusData[indBingo, c(1, 2, 4, 5)])
write.table(x = outCoor, file = outFile, quote = F, sep = "\t",
            row.names = F, col.names = T)


