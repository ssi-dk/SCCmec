#!/usr/bin/env python3

import os
import sys

#o = open(sys.argv[2],'w')
with open(sys.argv[1]) as f:
	for line in f:
		if line[0] == '>':
			splitline = line.rstrip('\n').split('[location')
			splitline = splitline[1].split(']')
			loc_string = splitline[0]
			ID = line.split(' ')[0]
			header = ID+'__loc'+loc_string
			print(header)
#o.close()