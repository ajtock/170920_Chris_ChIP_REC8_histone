# Determine peak midpoints and set start and end coordinates around midpoint
# so that all peaks have a common width of 200 bp
# Generate peak bed files (0-based start coordinates)

source("/projects/ajt200/Rfunctions/locus_midpoint.R")
library(GenomicRanges)

inDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/"
outDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/"

load(paste0(inDir, "armrangerPeaksGR_REC8_MYC_Rep1_minuslog10_p0.05_q0.05_qval_sorted_noMinWidth_1basedSummits.RData"))
seqlevels(armrangerPeaksGR) <- sub("Chr", "", seqlevels(armrangerPeaksGR))
armrangerPeaksGR_200bp <- locMidpointFlank(armrangerPeaksGR, leftFlank = 100, rightFlank = 100)
armrangerPeaks_200bp_0based <- cbind(as.numeric(seqnames(armrangerPeaksGR_200bp)),
                                     start(armrangerPeaksGR_200bp)-1, end(armrangerPeaksGR_200bp))
colnames(armrangerPeaks_200bp_0based) <- c("chr", "start", "end")
write.table(armrangerPeaks_200bp_0based, file = paste0(outDir, "armrangerPeaks_200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)

load(paste0(inDir, "perirangerPeaksGR_REC8_MYC_Rep1_minuslog10_p0.05_q0.05_qval_sorted_noMinWidth_1basedSummits.RData"))
seqlevels(perirangerPeaksGR) <- sub("Chr", "", seqlevels(perirangerPeaksGR))
perirangerPeaksGR_200bp <- locMidpointFlank(perirangerPeaksGR, leftFlank = 100, rightFlank = 100)
perirangerPeaks_200bp_0based <- cbind(as.numeric(seqnames(perirangerPeaksGR_200bp)),
                                      start(perirangerPeaksGR_200bp)-1, end(perirangerPeaksGR_200bp))
colnames(perirangerPeaks_200bp_0based) <- c("chr", "start", "end")
write.table(perirangerPeaks_200bp_0based, file = paste0(outDir, "perirangerPeaks_200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)


