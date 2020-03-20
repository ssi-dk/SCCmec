#!/usr/bin/env python3

import os
import sys

with open("ccrA_all_uniq.fasta") as f:
	for line in f:
		line = line.rstrip('\n')
		if line[0] == '>':
			print('>ccrA_'+line[1:])
		else:
			print(line)

with open("ccrB_all_uniq.fasta") as f:
	for line in f:
		line = line.rstrip('\n')
		if line[0] == '>':
			print('>ccrB_'+line[1:])
		else:
			print(line)

with open("ccrC_all_uniq.fasta") as f:
	for line in f:
		line = line.rstrip('\n')
		if line[0] == '>':
			print('>ccrC_'+line[1:])
		else:
			print(line)
