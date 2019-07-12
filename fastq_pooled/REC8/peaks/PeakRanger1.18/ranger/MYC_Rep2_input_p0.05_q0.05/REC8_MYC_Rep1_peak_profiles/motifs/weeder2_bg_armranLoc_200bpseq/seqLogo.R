# Convert position weight matrices (PWMs) to sequence logos using seqLogo Bioconductor package

library(seqLogo)

inDir <- "/home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_armranLoc_200bpseq/"

num_pwm <- as.numeric(system(paste0("ls -1 ", inDir, "MAT*.pwm | wc -l"), intern = T))

for(i in 1:num_pwm) {
  pwm <- read.table(file = system(paste0("ls ", inDir, "MAT", i, "_*.pwm"), intern = T), skip = 1, row.names = 1)
  print(head(pwm))
  postscript(file = paste0(inDir, "MAT", i, ".eps"))
  seqLogo(pwm)
  dev.off()
}

pwm_list <- list()
for(i in 1:num_pwm) {
  pwm_list[[i]] <- read.table(file = system(paste0("ls ", inDir, "MAT", i, "_*.pwm"), intern = T), skip = 1, row.names = 1)
}
pdf(file = paste0(inDir, "REC8_MYC_Rep1_peak_motifs_weeder2_bg_armranLoc_200bpseq.pdf"))
for(i in 1:num_pwm) {
  seqLogo(pwm_list[[i]])
}
dev.off()


