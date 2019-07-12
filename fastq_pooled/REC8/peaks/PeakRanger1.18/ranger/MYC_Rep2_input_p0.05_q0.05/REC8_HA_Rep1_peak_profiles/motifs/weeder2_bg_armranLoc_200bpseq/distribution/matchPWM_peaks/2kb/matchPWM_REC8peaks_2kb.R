# Plot distribution of motifs in 2-kb flanking regions around REC8 peak midpoints

#library(Biostrings)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(segmentSeq)
library(regioneR)
#library(zoo)
library(TTR)

peakDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/"
motifDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_HA_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/"
outDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_HA_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_peaks/2kb/"
plotDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_HA_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_peaks/2kb/plots/"

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

load(paste0(peakDir, "armrangerPeaksGR_REC8_HA_Rep1_minuslog10_p0.05_q0.05_qval_sorted_noMinWidth_1basedSummits.RData"))
seqlevels(armrangerPeaksGR) <- sub("Chr", "", seqlevels(armrangerPeaksGR))
peaksGR <- armrangerPeaksGR
ranLocGR <- randomizeRegions(peaksGR, genome = genome, mask = mask, per.chromosome = T, allow.overlaps = T)

chr1 <- Athaliana$Chr1
chr2 <- Athaliana$Chr2
chr3 <- Athaliana$Chr3
chr4 <- Athaliana$Chr4
chr5 <- Athaliana$Chr5
chr.list <- list()
chr.list[[1]] <- chr1
chr.list[[2]] <- chr2
chr.list[[3]] <- chr3
chr.list[[4]] <- chr4
chr.list[[5]] <- chr5



num_pwm <- as.numeric(system(paste0("ls -1 ", motifDir, "MAT*.pwm | wc -l"), intern = T))

for(i in 1:num_pwm) {
  pwm <- read.table(file = system(paste0("ls ", motifDir, "MAT", i, "_*.pwm"), intern = T), skip = 1, row.names = 1)
  pwm_transposed <- t(pwm)
  print(head(pwm))
  postscript(file = paste0(motifDir, "MAT", i, ".eps"))
  seqLogo(pwm)
  dev.off()
}

pwm_list <- list()
for(i in 1:num_pwm) {
  pwm_list[[i]] <- read.table(file = system(paste0("ls ", motifDir, "MAT", i, "_*.pwm"), intern = T), skip = 1, row.names = 1)
}


MOTIF_1_ABI3VP1_tnt.VRN1 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_1_ABI3VP1_tnt.VRN1_col_a_m1.txt", intern = T)
MOTIF_2_ABI3VP1_tnt.VRN1 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_2_ABI3VP1_tnt.VRN1_col_a_m1.txt", intern = T)
MOTIF_3_ABI3VP1_tnt.VRN1 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_3_ABI3VP1_tnt.VRN1_col_a_m1.txt", intern = T)
MOTIF_4_HMG_tnt.3XHMGBOX1 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_4_HMG_tnt.3XHMGBOX1_colamp_a_m1.txt", intern = T)
MOTIF_7_AP2EREBP_tnt.ERF115 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_7_AP2EREBP_tnt.ERF115_col_a_m1.txt", intern = T)
MOTIF_13_AP2EREBP_tnt.ERF48 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_13_AP2EREBP_tnt.ERF48_col_a_m1.txt", intern = T)
MOTIF_15_LOBAS2_tnt.LOB <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_15_LOBAS2_tnt.LOB_col_a_m1.txt", intern = T)
MOTIF_40_LOBAS2_tnt.LBD13 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_40_LOBAS2_tnt.LBD13_col_a_m1.txt", intern = T)
MOTIF_68_LOBAS2_tnt.LBD13 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_68_LOBAS2_tnt.LBD13_colamp_a_m1.txt", intern = T)
MOTIF_17_C2C2gata_tnt.ZML1 <- system("ls /projects/ajt200/REC8_MSH4/REC8_peaks/log2ChIPinput/nucleR/no_trim/analysis_03/motif_analysis/MEME-ChIP/MEME-ChIP_REC8armPeaksSH99merge_500bp_3o_m10/distribution/MOTIF_17_C2C2gata_tnt.ZML1_col_a_m1.txt", intern = T) 

mat.paths <- list(MOTIF_1_ABI3VP1_tnt.VRN1,
                  MOTIF_2_ABI3VP1_tnt.VRN1,
                  MOTIF_3_ABI3VP1_tnt.VRN1,
                  MOTIF_4_HMG_tnt.3XHMGBOX1,
                  MOTIF_7_AP2EREBP_tnt.ERF115,
                  MOTIF_13_AP2EREBP_tnt.ERF48,
                  MOTIF_15_LOBAS2_tnt.LOB,
                  MOTIF_40_LOBAS2_tnt.LBD13,
                  MOTIF_68_LOBAS2_tnt.LBD13,
                  MOTIF_17_C2C2gata_tnt.ZML1)
mat.names <- c("MOTIF_1_ABI3VP1_tnt.VRN1",
               "MOTIF_2_ABI3VP1_tnt.VRN1",
               "MOTIF_3_ABI3VP1_tnt.VRN1",
               "MOTIF_4_HMG_tnt.3XHMGBOX1",
               "MOTIF_7_AP2EREBP_tnt.ERF115",
               "MOTIF_13_AP2EREBP_tnt.ERF48",
               "MOTIF_15_LOBAS2_tnt.LOB",
               "MOTIF_40_LOBAS2_tnt.LBD13",
               "MOTIF_68_LOBAS2_tnt.LBD13",
               "MOTIF_17_C2C2gata_tnt.ZML1")
mat.list <- lapply(seq_along(mat.paths), function(x) {
  t(as.matrix(read.table(mat.paths[[x]])))
})
mat.list <- lapply(mat.list, function(x) {
  rownames(x) <- DNA_BASES; x
})
for(i in 1:length(mat.list)) {
  assign(paste0(mat.names[i]), mat.list[[i]])
}

mclapply(seq_along(mat.list), function(x) {
all.peaks.1201 <- rep(0, times = 1201)
all.ranLoc.1201 <- rep(0, times = 1201)
for(i in 1:5) {
  print(i)
  match.mat.list <- matchPWM(mat.list[[x]], chr.list[[i]], min.score = "80%")
  motif.ranges <- match.mat.list@ranges
  motif.starts <- motif.ranges@start
  # peaks
  print(i)
  chr.peaks <- peaksGR[seqnames(peaksGR) == i]
  mid.peaks <- chr.peaks$midpoint
  chr.peaks.1201 <- rep(0, times = 1201)  
  for(a in 1:length(mid.peaks)) {
    print(a)
    sel <- which(motif.starts >= mid.peaks[a]-600 & motif.starts <= mid.peaks[a]+600)
    motif.start.less.peak.start <- motif.starts[sel]-(mid.peaks[a]-601)
    sites <- rep(1, length(motif.start.less.peak.start))
    for(b in 1:length(motif.start.less.peak.start)) {
      chr.peaks.1201[motif.start.less.peak.start[b]] <- chr.peaks.1201[motif.start.less.peak.start[b]]+sites[b]
    }
  }
  all.peaks.1201 <- all.peaks.1201+chr.peaks.1201
  # ranLoc
  print(i)
  chr.ranLoc <- ranLocGR[seqnames(ranLocGR) == i]
  mid.ranLoc <- chr.ranLoc$midpoint 
  chr.ranLoc.1201 <- rep(0, times = 1201)
  for(a in 1:length(mid.ranLoc)) {
    print(a)
    sel <- which(motif.starts >= mid.ranLoc[a]-600 & motif.starts <= mid.ranLoc[a]+600)
    motif.start.less.ranLoc.start <- motif.starts[sel]-(mid.ranLoc[a]-601)
    sites <- rep(1, length(motif.start.less.ranLoc.start))
    for(b in 1:length(motif.start.less.ranLoc.start)) {
      chr.ranLoc.1201[motif.start.less.ranLoc.start[b]] <- chr.ranLoc.1201[motif.start.less.ranLoc.start[b]]+sites[b]
    }
  }
  all.ranLoc.1201 <- all.ranLoc.1201+chr.ranLoc.1201
  # Number of chr peaks
  print(paste0("Number of chromosome ", i, " peaks = ", length(mid.peaks)))
  # Number of chr ranLoc
  print(paste0("Number of chromosome ", i, " ranLoc = ", length(mid.ranLoc))) 
}
write.table(all.peaks.1201, file = paste0(outDir, mat.names[[x]], "_dist_REC8armPeaksSH99merged.txt"))
write.table(all.ranLoc.1201, file = paste0(outDir, mat.names[[x]],  "_dist_armRanLocMerged.txt"))
}, mc.cores = 10)

peaks.tab.list <- lapply(seq_along(mat.list), function(x) {
  read.table(file = paste0(outDir, mat.names[[x]], "_dist_REC8armPeaksSH99merged.txt"))
})
ranLoc.tab.list <- lapply(seq_along(mat.list), function(x) {
  read.table(file = paste0(outDir, mat.names[[x]], "_dist_armRanLocMerged.txt"))
})
peaks.tab.list.w20 <- lapply(seq_along(peaks.tab.list), function(x) {
  SMA(peaks.tab.list[[x]], 20)
})
ranLoc.tab.list.w20 <- lapply(seq_along(ranLoc.tab.list), function(x) {
  SMA(ranLoc.tab.list[[x]], 20)
})
peaks.tab.list.w30 <- lapply(seq_along(peaks.tab.list), function(x) {
  SMA(peaks.tab.list[[x]], 30)
})
ranLoc.tab.list.w30 <- lapply(seq_along(ranLoc.tab.list), function(x) {
  SMA(ranLoc.tab.list[[x]], 30)
})
peaks.tab.list.w40 <- lapply(seq_along(peaks.tab.list), function(x) {
  SMA(peaks.tab.list[[x]], 40)
})
ranLoc.tab.list.w40 <- lapply(seq_along(ranLoc.tab.list), function(x) {
  SMA(ranLoc.tab.list[[x]], 40)
})

motifDistPlot <- function(xplot, all.peaks.1201, all.ranLoc.1201, mycols, ylim, main) {
  # peaks
  plot(xplot, all.peaks.1201, col = mycols[1], lwd = 1.5, type = "l", ylim = ylim, ann = F, xaxt = "n")
  mtext(side = 2, line = 2, cex = 0.8, text = "Motif frequency")
  mtext(side = 3, line = 0.5, cex = 0.7, text = main)
  # ranLoc
  lines(xplot, all.ranLoc.1201, col = mycols[2], lwd = 1.5)
  axis(side = 1, at = c(-600, 0, 600), labels = c("-600 bp", "Midpoint", "+600 bp"))
  axis(side = 2, at = pretty(c(all.peaks.1201, all.ranLoc.1201)))
  abline(v = 0, lty = 3)
  box(lwd = 1.5)
}

pdf(paste0(plotDir, "motif_distributions_REC8armPeaksSH99merged_w20_v070917.pdf"), height = 12.5, width = 6)
par(mfrow = c(5,2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0))
#par(mar = c(5, 4, 4, 4))
#par(mgp = c(3, 0.75, 0))
xplot <- seq(-600, 600, by = 1)
mycols <- c("red", "black")
ylim <- lapply(seq_along(peaks.tab.list.w20), function(x) {
  c(min(c(peaks.tab.list.w20[[x]], ranLoc.tab.list.w20[[x]]), na.rm = T),
    max(c(peaks.tab.list.w20[[x]], ranLoc.tab.list.w20[[x]]), na.rm = T))
})
lapply(seq_along(peaks.tab.list.w20), function(x) {
  motifDistPlot(xplot, peaks.tab.list.w20[[x]], ranLoc.tab.list.w20[[x]], mycols,
                ylim = ylim[[x]], main = mat.names[x])
})
dev.off()

pdf(paste0(plotDir, "motif_distributions_REC8armPeaksSH99merged_w30_v070917.pdf"), height = 12.5, width = 6)
par(mfrow = c(5,2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0)) 
#par(mar = c(5, 4, 4, 4))
#par(mgp = c(3, 0.75, 0))
xplot <- seq(-600, 600, by = 1)
mycols <- c("red", "black")
ylim <- lapply(seq_along(peaks.tab.list.w30), function(x) {
  c(min(c(peaks.tab.list.w30[[x]], ranLoc.tab.list.w30[[x]]), na.rm = T), 
    max(c(peaks.tab.list.w30[[x]], ranLoc.tab.list.w30[[x]]), na.rm = T)) 
})
lapply(seq_along(peaks.tab.list.w30), function(x) {
  motifDistPlot(xplot, peaks.tab.list.w30[[x]], ranLoc.tab.list.w30[[x]], mycols,
                ylim = ylim[[x]], main = mat.names[x])
})
dev.off()

pdf(paste0(plotDir, "motif_distributions_REC8armPeaksSH99merged_w40_v070917.pdf"), height = 12.5, width = 6)
par(mfrow = c(5,2))
par(mar = c(2.1, 3.2, 2.1, 3.2))
par(mgp = c(2.25, 1, 0)) 
#par(mar = c(5, 4, 4, 4))
#par(mgp = c(3, 0.75, 0))
xplot <- seq(-600, 600, by = 1)
mycols <- c("red", "black")
ylim <- lapply(seq_along(peaks.tab.list.w40), function(x) {
  c(min(c(peaks.tab.list.w40[[x]], ranLoc.tab.list.w40[[x]]), na.rm = T), 
    max(c(peaks.tab.list.w40[[x]], ranLoc.tab.list.w40[[x]]), na.rm = T)) 
})
lapply(seq_along(peaks.tab.list.w40), function(x) {
  motifDistPlot(xplot, peaks.tab.list.w40[[x]], ranLoc.tab.list.w40[[x]], mycols,
                ylim = ylim[[x]], main = mat.names[x])
})
dev.off()

