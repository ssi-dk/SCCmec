#!/usr/bin/env python3

import os
import sys

IDs = []
with open("Meningococci_ID_list.txt") as f:
	for line in f:
		IDs.append(line.rstrip('\n'))

no_match = []
found = []
missing = IDs

with open("read_list.txt") as f:
	for line in f:
		ID = line.rstrip('\n')[:-21]
		if ID in IDs:
			found.append(ID)
			missing.remove(ID)
		else:
			no_match.append(ID)


print("\n".join(no_match))
