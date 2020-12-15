#!/bin/python
#Alan E. Yocca
#nucleotide_diversity_calc.py
#A little clunky since on time crunch, just looping same instances used for pin/pis calcs

import dendropy

#loop files
#load into array
file_list = []

wkdir="/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/32_diversity/pin_pis/"
with open(wkdir + "/edit_aln_list.txt") as fh:
	for line in fh:
		file_list.append(line.strip())

output = []
#gene	class	value
iter = 0
for file in file_list:
	iter+=1
	if iter % 1000 == 0:
		print("Through %s alignments" % (iter))
	
	seqs = dendropy.DnaCharacterMatrix.get(path=file,schema = "fasta")
	div = dendropy.calculate.popgenstat.nucleotide_diversity(seqs)
	gene = file.split("/")[-1].split("_")[0]
	cns_class = file.split("/")[-2].split("_")[1]
	output.append([gene,cns_class,div])

#/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/32_diversity/pin_pis/aln_list.txt
#36271 /mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/32_diversity/pin_pis/aln_list.txt
#so thats all instances with at least 2 in each class
#/mnt/research/edgerpat_lab/AlanY/01_athal_cns/assembly_golicz16/32_diversity/pin_pis/01_split_fasta/01_same/AT1G01010_aln.fasta
#Can't seem to find the ideal ones, but can just filter those later in R right? yea so just do these

with open("div_calc_by_class.txt", "w") as out:
	for line in output:
		iter = 1
		for item in line:
			if iter != 1:
				tmp = out.write("\t")
			iter = 0				
			tmp = out.write(item)
		tmp = out.write("\n")
