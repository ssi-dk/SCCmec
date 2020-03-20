#!/usr/bin/env python3

import os
import sys

def reverse_complement(in_seq):
	out_seq = in_seq[::-1].lower()
	out_seq = out_seq.replace('a','b')
	out_seq = out_seq.replace('t','a')
	out_seq = out_seq.replace('b','t')
	out_seq = out_seq.replace('c','b')
	out_seq = out_seq.replace('g','c')
	out_seq = out_seq.replace('b','g')
	out_seq = out_seq.upper()
	return out_seq;

blast_file = sys.argv[1]
assembly_dir = sys.argv[2]
cov_req = float(sys.argv[3])
pident_req = float(sys.argv[4])
start_gap = int(sys.argv[5])
end_gap = int(sys.argv[6])

assembly_dir = '/srv/data/DB/refseq/Staphylococcus_combined'

already_tested = os.listdir('/srv/data/MPV/THEJ/Projekter/spl/spl_operon_combined_prokka/')
already_tested = []

blast_dict = []

with open(blast_file) as f:
	for line in f:
		line = line.rstrip('\n').split('\t')
		qlen = int(line[0])
		qstart = int(line[1])
		qend = int(line[2])
		pident = float(line[3])
		ID = line[4]
		sstart = int(line[7])
		send = int(line[8])
		if (qend-qstart+1)>=(qlen*cov_req/100) and pident>=pident_req:
			blast_dict.append([ID,qstart,qend,sstart,send])
			
for l in blast_dict:
	ID = l[0]
	qs = l[1]
	qe = l[2]
	ss = l[3]
	se = l[4]
	if ss < se:
		contig_start = int(ss-qs-start_gap)
		contig_end = int(se+qlen-qe+end_gap)
		rev_flag = 0
	else:
		contig_start = int(se-qs-end_gap)
		contig_end = int(ss+qlen-qe+start_gap)
		rev_flag = 1
	IDsplit = ID.split('|')
	assembly = IDsplit[0]
	contig = IDsplit[1]
	path = os.path.join(assembly_dir,assembly+'_genomic.fna')
	check_name = assembly+'__'+contig+'.fasta'
	if check_name in already_tested:
		print("Skipping "+check_name, file=sys.stderr)
	else:
		if os.path.exists:
			print("Checking "+check_name, file=sys.stderr)
			with open(path) as f:
				flag = 0
				for line in f:
					line = line.rstrip('\n')
					if line[0] == '>':
						if flag == 1:
							#print('\t'.join([ID,str(contig_start),str(contig_end),str(rev_flag)]))
							if contig_start < 1:
								contig_start = 1
							if contig_end > len(contig_seq):
								contig_end = len(contig_seq)
							if rev_flag == 0:
								seq = contig_seq[contig_start:contig_end]
							else:
								seq = reverse_complement(contig_seq[contig_start:contig_end])
							print('>'+ID+'\n'+seq)
						flag = 0
						header = line.split(' ')
						if header[0][1:] == contig:
							flag = 1
							contig_seq = ''
					else:
						if flag == 1:
							contig_seq += line
				if flag == 1:
					#print('\t'.join([ID,str(contig_start),str(contig_end),str(rev_flag)]))
					if contig_start < 1:
						contig_start = 1
					if contig_end > len(contig_seq):
						contig_end = len(contig_seq)
					if rev_flag == 0:
						seq = contig_seq[contig_start:contig_end]
					else:
						seq = reverse_complement(contig_seq[contig_start:contig_end])
					print('>'+ID+'\n'+seq)

					