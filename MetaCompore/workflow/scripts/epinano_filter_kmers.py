#!/bin/python

import sys

r = ["A", "G"]
h = ["A", "C", "T", "U"]

input = snakemake.input.variants
output= snakemake.output.filteredvariants

with open(output,'w') as outfile:
    with open(input, 'r') as f:
        for line in f:
            line=line
            if line[0]=="#":
                outfile.write(line)
            else:
                fields=line.split(',')
                kmer=fields[0].upper()
                if kmer[0] in r and kmer[1] in r and kmer[2] == "A" and kmer[3]=="C" and kmer[4] in h:
                    outfile.write(line)
            
