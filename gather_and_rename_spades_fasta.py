#!/usr/bin/env python3

import os
import sys

in_dir = sys.argv[1]
out_dir = sys.argv[2]

files  = os.listdir(in_dir)

if not os.path.exists(out_dir):
	os.makedirs(out_dir)

for file in files:
	file_split = file.split('_')
	ID = '_'.join(file_split[:-2])
	in_file = os.path.join(in_dir,file,'contigs.fasta')
	printline = ''
	if os.path.exists(in_file):
		with open(in_file) as f:
			for line in f:
				if line[0] == '>':
					printline += '>'+ID+'__'+line[1:]
				else:
					printline += line
		out_file = os.path.join(out_dir,ID+'.fasta')
		o = open(out_file,'w')
		o.write(printline)
		o.close()
	else:
		print(in_file)
