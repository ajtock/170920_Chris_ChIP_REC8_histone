# Determine peak midpoints and set start and end coordinates around midpoint
# so that all peaks have a common width of 200 bp
# Define random loci based on these 200-bp peak loci
# Generate random loci bed files (0-based start coordinates)

source("/projects/ajt200/Rfunctions/locus_midpoint.R")
library(GenomicRanges)
library(regioneR)

chrs <- c("Chr1","Chr2","Chr3","Chr4","Chr5")
chrStart <- c(1, 1, 1, 1, 1)
chrLens <- c(30427671, 19698289, 23459830, 18585056, 26975502)
centromeres <- c(15086045, 3607929, 13587786, 3956021, 11725024)
pericenStart <- c(11330001, 990001, 10200001, 990001, 8890001)
pericenEnd <- c(18480000, 7540000, 16860000, 6850000, 15650000)
genome <- toGRanges(data.frame(chrs, chrStart, chrLens))
seqlevels(genome) <- sub("Chr", "", seqlevels(genome))

inDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/"
outDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_HA_Rep1_peak_profiles/motifs/"

## arm
load(paste0(inDir, "armrangerPeaksGR_REC8_HA_Rep1_minuslog10_p0.05_q0.05_qval_sorted_noMinWidth_1basedSummits.RData"))
seqlevels(armrangerPeaksGR) <- sub("Chr", "", seqlevels(armrangerPeaksGR))
armrangerPeaksGR_200bp <- locMidpointFlank(armrangerPeaksGR, leftFlank = 100, rightFlank = 100)
armrangerPeaksGR_200bp_fmt <- GRanges(seqnames = seqnames(armrangerPeaksGR_200bp), ranges = ranges(armrangerPeaksGR_200bp), strand = strand(armrangerPeaksGR_200bp))

# mask pericentromeric regions
maskPeri <- toGRanges(data.frame(chrs, pericenStart, pericenEnd))
seqlevels(maskPeri) <- sub("Chr", "", seqlevels(maskPeri))
armranLocGR_200bp <- randomizeRegions(armrangerPeaksGR_200bp, genome = genome, mask = c(maskPeri, armrangerPeaksGR_200bp_fmt), per.chromosome = TRUE, allow.overlaps = TRUE)

armranLoc_200bp_0based <- cbind(as.numeric(seqnames(armranLocGR_200bp)),
                                start(armranLocGR_200bp)-1, end(armranLocGR_200bp))
colnames(armranLoc_200bp_0based) <- c("chr", "start", "end")
write.table(armranLoc_200bp_0based, file = paste0(outDir, "armranLoc_200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)


## peri
load(paste0(inDir, "perirangerPeaksGR_REC8_HA_Rep1_minuslog10_p0.05_q0.05_qval_sorted_noMinWidth_1basedSummits.RData"))
seqlevels(perirangerPeaksGR) <- sub("Chr", "", seqlevels(perirangerPeaksGR))
perirangerPeaksGR_200bp <- locMidpointFlank(perirangerPeaksGR, leftFlank = 100, rightFlank = 100)
perirangerPeaksGR_200bp_fmt <- GRanges(seqnames = seqnames(perirangerPeaksGR_200bp), ranges = ranges(perirangerPeaksGR_200bp), strand = strand(perirangerPeaksGR_200bp))

# mask arm regions
maskArms <- toGRanges(data.frame(rep(chrs, 2), c(chrStart, pericenEnd), c(pericenStart, chrLens)))
seqlevels(maskArms) <- sub("Chr", "", seqlevels(maskArms))
periranLocGR_200bp <- randomizeRegions(perirangerPeaksGR_200bp, genome = genome, mask = c(maskArms, perirangerPeaksGR_200bp_fmt), per.chromosome = TRUE, allow.overlaps = TRUE)

periranLoc_200bp_0based <- cbind(as.numeric(seqnames(periranLocGR_200bp)),
                                 start(periranLocGR_200bp)-1, end(periranLocGR_200bp))
colnames(periranLoc_200bp_0based) <- c("chr", "start", "end")
write.table(periranLoc_200bp_0based, file = paste0(outDir, "periranLoc_200bp_0based.bed"),
            row.names = F, col.names = F, quote = F)


