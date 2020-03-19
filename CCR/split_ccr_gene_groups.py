#!/usr/bin/env python3

import os
import sys


g_dict = {}
with open(sys.argv[1]) as f:
	for line in f:
		line = line.rstrip('\n').split('\t')
		g_dict[line[0]] = [line[8],line[11]]

fasta_dict = []
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
		print(header)