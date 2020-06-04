#!/usr/bin/env python3

import os
import sys

in_dir = "/srv/data/MPV/MTG/meningococci_Orebro/reads"
out_dir = "/srv/data/MPV/MTG/meningococci_Orebro/reads_renamed"
rename_list = "meningococci_rename.txt"
rename_dict = {}

with open(rename_list) as f:
	for line in f:
		line = line.rstrip('\n').split('\t')
		rename_dict[line[1]] = line[0]

files = os.listdir(in_dir)

for file in files:
	s = file.split('_')
	ID = '_'.join(s[:-3])
	suffix = '_'.join(s[-2:])
	new_ID = rename_dict[ID]
	src = os.path.join(in_dir,file)
	dst = os.path.join(out_dir,new_ID+'_'+suffix)
	os.symlink(src,dst)