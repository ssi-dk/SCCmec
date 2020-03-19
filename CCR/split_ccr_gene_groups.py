#!/usr/bin/env python3

import os
import sys


g_dict = {}
with open(sys.argv[1]) as f:
	for line in f:
		line = line.rstrip('\n').split('\t')
		g_dict[line[0]] = [line[8],line[11]]

fasta_dict = {}
with open(sys.argv[2]) as f:
	for line in f:
		line = line.rstrip('\n')
		if line[0] == '>':
			header = line[1:]
			fasta_dict[header] = ''
		else:
			fasta_dict[header] += line

print_dict = {}
for header in fasta_dict:
	if header in g_dict:
		allotype = g_dict[header][1]
		if allotype in print_dict:
			print_dict[allotype] += '>'+header+'\n'+fasta_dict[header]+'\n'
		else:
			print_dict[allotype] = '>'+header+'\n'+fasta_dict[header]+'\n'

for allotype in print_dict:
	out_file = os.path.join(sys.argv[3],allotype+'.fasta')
	o = open(out_file,'w')
	o.write(print_dict[allotype])
	o.close()