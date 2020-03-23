#!/usr/bin/env python3
import os
import sys


ignore_file = sys.argv[1]
fasta_file = sys.argv[2]

ignore_list = []
with open(ignore_file) as f:
	for line in f:
		line = line.rstrip('\n').split('__')[1]
		ignore_list.append(line)


with open(fasta_file) as f:
	for line in f:
		line = line.rstrip('\n')
		if line[0] == '>':
			test = line.split('__')[1]
			if test in ignore_list:
				flag = 0
			else:
				flag = 1
		if flag == 1:
			print(line)