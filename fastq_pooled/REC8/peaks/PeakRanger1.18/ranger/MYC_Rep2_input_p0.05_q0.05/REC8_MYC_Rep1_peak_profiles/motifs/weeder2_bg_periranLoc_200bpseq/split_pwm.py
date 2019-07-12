#!/usr/bin/env python

# GENERAL USAGE: python ./split_pwm.py ./perirangerPeaks_200bpseq.fa.matrix.w2

import os
from os import path
import sys

infile = open(sys.argv[1])

os.system("cd /home/meiosis/ajt200/analysis/170920_Chris_ChIP_REC8_histone/fastq_pooled/REC8/peaks/PeakRanger1.18/ranger/MYC_Rep2_input_p0.05_q0.05/REC8_MYC_Rep1_peak_profiles/motifs/weeder2_bg_periranLoc_200bpseq") # Add project directory here

path = sys.argv[1].split('/')[0]

matrix_file = path + "/" + sys.argv[1]
print('\n' + "matrix file: " + matrix_file + '\n')
num_pwms = 'grep "^>" %s | wc -l' % (matrix_file)
print("Number of pwms in matrix file = ")
os.system(num_pwms)

opened = False # Assume outfile is not open

for line_matrix in infile:
    if line_matrix[0] == ">": # If line begins with ">"
        if(opened):
            outfile.close() # Will close the outfile if it is open (see below and follow loop)
        opened = True # Set opened to True to represent an opened outfile
        pwm_name = line_matrix[1:].rstrip().replace("\t", "_") # Extract pwm name: remove ">", extract pwm string, remove any spaces or new lines following file
        #pwm_name = line_matrix[1:6].rstrip()
        print("pwm: " + pwm_name)
        print("Name of outfile: " + path + "/" + str(pwm_name) + ".pwm")
        outfile = open(path + "/" + str(pwm_name) + ".pwm", 'w')
    #if line_matrix[0] != ">":
    outfile.write(line_matrix) # Write the line to the file. If the line starts with ">", a new output file will be created, ready to be written to.
                               # Subsequent lines (i.e., the base relative frequencies) will all be written as the script loops through lines.
                               # Only when the script reaches a new line beginning with ">" is the outfile closed,
                               # and a new outfile with a new name is generated.
outfile.close()
print "Finished"

