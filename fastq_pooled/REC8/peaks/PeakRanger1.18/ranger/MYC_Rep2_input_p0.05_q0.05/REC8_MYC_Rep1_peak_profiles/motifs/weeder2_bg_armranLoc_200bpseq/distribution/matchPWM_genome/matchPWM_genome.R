# Identify genomic loci that match motifs enriched at REC8 peaks
# Profile REC8 and nucleosomes at these loci
# Profile base composition (proportion and -log10(probability)) in regions flanking these loci

library(Biostrings)
library(BSgenome.Athaliana.TAIR.TAIR9)
library(segmentSeq)
library(regioneR)
#library(zoo)
library(TTR)

motifDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/"
outDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/"
plotDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/distribution/matchPWM_genome/plots/"

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

chr1 <- Athaliana$Chr1
chr2 <- Athaliana$Chr2
chr3 <- Athaliana$Chr3
chr4 <- Athaliana$Chr4
chr5 <- Athaliana$Chr5
chr_list <- list()
chr_list[[1]] <- chr1
chr_list[[2]] <- chr2
chr_list[[3]] <- chr3
chr_list[[4]] <- chr4
chr_list[[5]] <- chr5

num_pwm <- as.numeric(system(paste0("ls -1 ", motifDir, "MAT*.pwm | wc -l"), intern = T))

pwm_list <- list()
for(i in 1:num_pwm) {
  pwm_list[[i]] <- as.matrix(read.table(file = system(paste0("ls ", motifDir, "MAT", i, "_*.pwm"), intern = T), skip = 1, row.names = 1))
}

mclapply(seq_along(pwm_list), function(x) {
  motif.GRanges <- GRanges()
  for(i in 1:5) {
    print(i)
    match.pwm_list <- matchPWM(pwm_list[[x]], chr_list[[i]], min.score = "87.5%")
    motif.ranges <- match.pwm_list@ranges
    motif.GRanges.chr <- GRanges(seqnames = i, ranges = motif.ranges, strand = "*")
    motif.GRanges <- append(motif.GRanges, motif.GRanges.chr)
    }
    save(motif.GRanges, file = paste0(outDir, "motif", x, "_matchPWM_GRanges.RData"))
}, mc.cores = length(pwm_list))


library(doParallel)
# Change number of cores to reflect number of samples you want to process simultaneously
registerDoParallel(cores = length(pwm_list))
print("Currently registered parallel backend name, version and cores")
print(getDoParName())
print(getDoParVersion())
print(getDoParWorkers())

foreach(x = 1:num_pwm) %dopar% {
  print(x)
  load(file = paste0(outDir, "motif", x, "_matchPWM_GRanges.RData"))
  tmp <- DNAStringSet()
  for(h in 1:5) {
    tmp.chr <- DNAStringSet()
    # Obtain sequences for each motif match and flanking 20-bp regions as DNAStringSet object
    for(i in 1:length(motif.GRanges[seqnames(motif.GRanges) == h])) {
      tmp.chr <- c(tmp.chr, DNAStringSet(chr_list[[h]][(start(ranges(motif.GRanges)[seqnames(motif.GRanges) == h][i])-20):(end(ranges(motif.GRanges)[seqnames(motif.GRanges) == h][i])+20)]))
    }
    tmp <- c(tmp, tmp.chr)
  }
  # Generate position frequency matrix (PFM)
  pfm <- consensusMatrix(tmp)
  # Convert frequencies to proportions and retain rows 1:4
  prm <- prop.table(pfm, 2)[1:4,]
  # Re-order rows for stack barplot representation
  prm_AGTC <- rbind(prm[1,], prm[3,], prm[4,], prm[2,])
  rownames(prm_AGTC) <- c("A", "G", "T", "C")
  pdf(paste0(plotDir, "REC8_MYC_Rep1_peak_motif", x, "_", paste0(strsplit(consensusString(tmp), split = "")[[1]][21:(21+mean(width(motif.GRanges))-1)], collapse = ""), "_base_proportions.pdf"))
  par(mgp = c(2, 1, 0))
  barplot(prm_AGTC,
          col = c("green", "yellow", "red", "blue"),
          xlab = paste0("Position within REC8-MYC Rep1 peak motif", x, "_", paste0(strsplit(consensusString(tmp), split = "")[[1]][21:(21+mean(width(motif.GRanges))-1)], collapse = ""), " matches and 20-bp flanks"),
          ylab = "Proportion",
          legend.text = TRUE,
          args.legend = list(
          x = ncol(prm_AGTC) + 16,
          y = 0.2,
          bty = "n"
          )
         )
  dev.off()
}

