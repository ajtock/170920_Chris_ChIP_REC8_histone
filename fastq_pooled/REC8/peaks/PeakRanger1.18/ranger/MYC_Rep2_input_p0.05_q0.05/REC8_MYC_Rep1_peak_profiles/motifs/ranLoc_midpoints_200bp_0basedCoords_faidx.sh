#!/bin/bash

# obtain locus sequences in fasta format with coordinates in 0-based bed file

faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed armranLoc_200bp_0based.bed --out armranLoc_200bpseq.fa
faidx /projects/ajt200/TAIR10/TAIR10_chr_all.fa --bed periranLoc_200bp_0based.bed --out periranLoc_200bpseq.fa

