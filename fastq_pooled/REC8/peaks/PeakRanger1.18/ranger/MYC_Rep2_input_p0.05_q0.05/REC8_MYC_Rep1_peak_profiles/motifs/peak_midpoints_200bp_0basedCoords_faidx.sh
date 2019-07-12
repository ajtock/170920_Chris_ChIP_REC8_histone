#!/bin/bash

# obtain locus sequences in fasta format with coordinates in 0-based bed file

faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed armrangerPeaks_200bp_0based.bed --out armrangerPeaks_200bpseq.fa
faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed perirangerPeaks_200bp_0based.bed --out perirangerPeaks_200bpseq.fa

