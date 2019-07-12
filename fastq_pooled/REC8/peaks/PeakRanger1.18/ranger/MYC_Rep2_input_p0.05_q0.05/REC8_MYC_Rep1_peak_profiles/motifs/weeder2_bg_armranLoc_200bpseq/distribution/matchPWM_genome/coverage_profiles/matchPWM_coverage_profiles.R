# Profile REC8 and nucleosomes at genomic loci that match motifs enriched at REC8 peaks

library(segmentSeq)
library(parallel)
library(EnrichedHeatmap)
library(genomation)
library(regioneR)

motifDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/"
lociDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/"
matDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/coverage_profiles/matrices/"
plotDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/coverage_profiles/plots/"

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
seqlevels(genome) <- sub("Chr", "", seqlevels(genome))
mask <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))
seqlevels(mask) <- sub("Chr", "", seqlevels(mask))

# specify locations of normalised per base coverage files
REC8_HA_Rep1 <- system("ls /home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/coverage/common_input_MYC_Rep2/log2ChIPinput/noZscore/log2_REC8_HA_Rep1_ChIP_MYC_Rep2_input_noZscore_norm_allchrs_coverage_coord_tab.bed", intern = T)
REC8_MYC_Rep2 <- system("ls /home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/coverage/common_input_MYC_Rep2/log2ChIPinput/noZscore/log2_REC8_MYC_Rep2_ChIP_MYC_Rep2_input_noZscore_norm_allchrs_coverage_coord_tab.bed", intern = T)
REC8_MYC_Rep1 <- system("ls /home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/coverage/common_input_MYC_Rep2/log2ChIPinput/noZscore/log2_REC8_MYC_Rep1_ChIP_MYC_Rep2_input_noZscore_norm_allchrs_coverage_coord_tab.bed", intern = T)
MNase <- system("ls /projects/ajt200/BAM_masters/nucleosomes/WT/coverage/nakedDNA_untrimmed_input/log2ChIPinput/noZscore/log2wtNucNakedDNAuntrimmed_noZscore_norm_allchrs_coverage_coord_tab.bed", intern = T)

# create lists/vectors of path and library names
libPaths <- list(REC8_HA_Rep1, REC8_MYC_Rep2, REC8_MYC_Rep1, MNase)
libNames <- c("REC8_HA_Rep1", "REC8_MYC_Rep2", "REC8_MYC_Rep1", "MNase")

# import coverage files as GRanges objects and assign to library names

# An unfortunate feature of readGeneric() is that if the meta.col (coverage column) begins with a series of "0"s,
# the column is automatically determined to be of class "integer" rather than of the desired class "numeric".
# So it is necessary to change (e.g., with vim) the first coverage value from "0" to "0.0" in the problem BED-like file.
grTmp <- mclapply(seq_along(libPaths), function(x) {
  #import(libPaths[[x]], format = "bed")
  readGeneric(libPaths[[x]], meta.col = list(coverage = 4))
}, mc.cores = 4, mc.preschedule = F)
for(i in 1:length(grTmp)) {
  seqlevels(grTmp[[i]]) <- sub("Chr", "", seqlevels(grTmp[[i]]))
  #grTmp[[i]] <- GRanges(seqnames = seqnames(grTmp[[i]]), ranges = IRanges(start = start(grTmp[[i]])-1, end = end(grTmp[[i]])),
  #                      strand = strand(grTmp[[i]]), coverage = as.numeric(grTmp[[i]]$name))
  assign(paste0(libNames[i]), grTmp[[i]])
}

# Create GRangesList object containing per base coverage for each library
grl <- GRangesList("REC8_HA_Rep1" = REC8_HA_Rep1, "REC8_MYC_Rep2" = REC8_MYC_Rep2, "REC8_MYC_Rep1" = REC8_MYC_Rep1, "MNase" = MNase)

# Function to create coverage matrices for target loci and random loci (incl. flanking regions)
## and to calculate mean levels per window
covMatrix <- function(signal, target, ranLoc, y) {
  #target loci
  set.seed(2840)
  mat1 <- normalizeToMatrix(signal, target, value_column = "coverage",
                            extend = flankSize, mean_mode = "absolute", w = winSize,
                            empty_value = 0, smooth = FALSE,
                            include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat1)
  print(length(mat1))
  mat1_DF <- data.frame(mat1)
  mat1_DF_colMeans <- as.vector(colMeans(mat1_DF))
  write.table(mat1_DF_colMeans, file = outDFCM[[y]][[1]])

  #random loci
  set.seed(8472)
  mat2 <- normalizeToMatrix(signal, ranLoc, value_column = "coverage",
                            extend = flankSize, mean_mode = "absolute", w = winSize,
                            empty_value = 0, smooth = FALSE,
                            include_target = TRUE, target_ratio = targetSize/(targetSize+(flankSize*2)))
  print(mat2)
  print(length(mat2))
  mat2_DF <- data.frame(mat2)
  mat2_DF_colMeans <- as.vector(colMeans(mat2_DF))
  write.table(mat2_DF_colMeans, file = outDFCM[[y]][[2]])
}

# Function to plot mean coverage profiles of each REC8 replicate vs MNase around motif matches and random loci
plotCov_motifMatch <- function(xplot, dat1, dat2, ranDat1, ranDat2, Ylabel1, Ylabel2, x) {
  plot(xplot, dat1[,1],
       ylim = c(min(dat1[,1], ranDat1[,1]),
                max(dat1[,1], ranDat1[,1])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 1, cex = 0.8, text = paste("REC8-MYC Rep1 arm peak motif", x, " matched loci", sep = ""))
  axis(side = 2, at = pretty(c(dat1[,1], ranDat1[,1])))
  mtext(side = 2, line = 2, cex = 0.8, text = Ylabel1, col = mycols[1])
  par(new = T)
  plot(xplot, dat2[,1], ylim = c(min(dat2[,1], ranDat2[,1]),
                                 max(dat2[,1], ranDat2[,1])),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2[,1], ranDat2[,1])))
  axis(side = 1, at = c(0, flankSize, length(dat1[,1])-flankSize, length(dat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(0, flankSize, length(dat1[,1])-flankSize, length(dat1[,1])), text = c("-1 kb", "Start", "End", "+1 kb"))
  abline(v = c(flankSize, length(dat1[,1])-flankSize), lty = 3)
  box(lwd = 1.5)

  plot(xplot, ranDat1[,1], 
       ylim = c(min(dat1[,1], ranDat1[,1]),
                max(dat1[,1], ranDat1[,1])),
       type = "l", lwd = 1.5, col = mycols[1], ann = F, xaxt = "n", yaxt = "n")
  mtext(side = 3, line = 1, cex = 0.8, text = "Random loci")
  axis(side = 2, at = pretty(c(dat1[,1], ranDat1[,1])))
  par(new = T)
  plot(xplot, ranDat2[,1],
       ylim = c(min(dat2[,1], ranDat2[,1]),
                max(dat2[,1], ranDat2[,1])),
       type = "l", lwd = 1.5, col = mycols[2], ann = F, xaxt = "n", yaxt = "n")
  axis(side = 4, at = pretty(c(dat2[,1], ranDat2[,1])))
  mtext(side = 4, line = 2, cex = 0.8, text = Ylabel2, col = mycols[2])
  axis(side = 1, at = c(0, flankSize, length(ranDat1[,1])-flankSize, length(ranDat1[,1])), labels = c("", "", "", ""))
  mtext(side = 1, line = 1, cex = 0.7, at = c(0, flankSize, length(ranDat1[,1])-flankSize, length(ranDat1[,1])), text = c("-1 kb", "Start", "End", "+1 kb"))
  abline(v = c(flankSize, length(ranDat1[,1])-flankSize), lty = 3)
  box(lwd = 1.5)
}

num_pwm <- as.numeric(system(paste0("ls -1 ", motifDir, "MAT*.pwm | wc -l"), intern = T))

#library(doParallel)
## Change number of cores to reflect number of samples you want to process simultaneously
#registerDoParallel(cores = num_pwm)
#print("Currently registered parallel backend name, version and cores")
#print(getDoParName())
#print(getDoParVersion())
#print(getDoParWorkers())

for(x in 1:num_pwm) {
  print(x)
  load(file = paste0(lociDir, "motif", x, "_matchPWM_GRanges.RData"))
  motifMatchGR <- motif.GRanges
  ranLocGR <- randomizeRegions(motifMatchGR, genome = genome, mask = mask, per.chromosome = TRUE, allow.overlaps = TRUE)

  winSize <- 1
  targetSize <- mean(width(motifMatchGR))
  flankSize <- 1000

  # Define column mean outfiles
  outDFCM <- lapply(seq_along(libNames), function(y)
    list(paste0(matDir, libNames[[y]],
                "_norm_cov_motif", x, "_matches_mat1_target_and_flank_dataframe_colMeans.txt"),
         paste0(matDir, libNames[[y]],
                "_norm_cov_ranLoc", x, "_mat2_target_and_flank_dataframe_colMeans.txt")))

  # Run covMatrix function on each coverage dataset
  mclapply(seq_along(grl), function(y) {
    covMatrix(grl[[y]], motifMatchGR, ranLocGR, y)
  }, mc.cores = 4, mc.preschedule = F)

  # Load column mean outfiles
  motifMatchGR_dat1 <- read.table(file = outDFCM[[1]][[1]])
  motifMatchGR_dat2 <- read.table(file = outDFCM[[2]][[1]])
  motifMatchGR_dat3 <- read.table(file = outDFCM[[3]][[1]])
  motifMatchGR_dat4 <- read.table(file = outDFCM[[4]][[1]])
  ranLocGR_dat1 <- read.table(file = outDFCM[[1]][[2]])
  ranLocGR_dat2 <- read.table(file = outDFCM[[2]][[2]])
  ranLocGR_dat3 <- read.table(file = outDFCM[[3]][[2]])
  ranLocGR_dat4 <- read.table(file = outDFCM[[4]][[2]])

  # Plot mean REC8 and nucleosomes coverage profiles around genomic loci that match motifs enriched at REC8 peaks
  # and around random loci
  pdf(paste0(plotDir, "REC8_MYC_Rep1_peak_motif", x, "_matched_loci_coverage_profiles_REC8_MNase.pdf"), height = 7.5, width = 6)
  par(mfrow = c(3, 2))
  par(mar = c(2.1, 3.2, 2.1, 3.2))
  par(mgp = c(2.25, 1, 0))
  xplot <- seq(1, length(motifMatchGR_dat1[,1]), by = 1)
  mycols <- c("blue", "red")

  plotCov_motifMatch(xplot = xplot, dat1 = motifMatchGR_dat1, dat2 = motifMatchGR_dat4,
                     ranDat1 = ranLocGR_dat1, ranDat2 = ranLocGR_dat4,
                     Ylabel1 = "REC8-HA Rep1", Ylabel2 = "MNase", x)
  plotCov_motifMatch(xplot = xplot, dat1 = motifMatchGR_dat2, dat2 = motifMatchGR_dat4,
                     ranDat1 = ranLocGR_dat2, ranDat2 = ranLocGR_dat4,
                     Ylabel1 = "REC8-MYC Rep2", Ylabel2 = "MNase", x)
  plotCov_motifMatch(xplot = xplot, dat1 = motifMatchGR_dat3, dat2 = motifMatchGR_dat4,
                     ranDat1 = ranLocGR_dat3, ranDat2 = ranLocGR_dat4,
                     Ylabel1 = "REC8-MYC Rep1", Ylabel2 = "MNase", x)
  dev.off()
}


