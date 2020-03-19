#!/usr/bin/env python3

import os
import sys
import numpy as np

fasta_dict = {}
headers = []
with open(sys.argv[1]) as f:
	for line in f:
		line = line.rstrip('\n')
		if line[0] == '>':
			header = line[1:]
			fasta_dict[header] = ''
			headers.append(header)
		else:
			fasta_dict[header] += line

sample_seq = fasta_dict[header]

dist_matrix = np.zeros((len(headers),len(headers)))


for i in range(len(headers)):
	seq_1 = fasta_dict[headers[i]]
	seq_1_len = len(seq_1.replace('-','')
	for j in range(len(headers)):
		diff_counts = 0
		seq_2= fasta_dict[headers[j]]
		seq_2_len = len(seq_2.replace('-','')
		for n in range(len(seq_1)):
			if seq_1[n] != seq_2[n] and seq_1[n].upper() != 'N' and seq_2[n].upper() != 'N' and seq_1[n] != '-' and seq_2[n] != '-':
				diff_counts += 1
		seq_len = min(seq_1_len,seq_2_len)
		dist_matrix[i][j] = diff_counts/seq_len*100
		dist_matrix[j][i] = diff_counts/seq_len*100



print('\t'+'\t'.join(headers))
for i in range(len(headers)):
	printlist = [headers[i]]
	for j in range(len(headers)):
		printlist.append(str(int(dist_matrix[i][j])))
	print('\t'.join(printlist))
	