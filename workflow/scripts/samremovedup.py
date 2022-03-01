#!/usr/bin/env python3

"""
Author: Pontus Skoglund
pontus.skoglund@gmail.com

Script modifications: Verena Kutschera
"""

import sys

lastRead = None
lastlengths=[]
for line in sys.stdin:
	if line[0] == '@':
		print(line, end="")
		continue
		
	fields = line.split()
	flag = int(fields[1])
	if flag & 0x4: #unmapped read
		print(line, end="")
		continue #Verena: print also unmapped reads to output
    
	chrname = fields[2]
	pos = int(fields[3])
	seq = fields[9]
	length=len(seq)
	
	if lastRead == (chrname, pos) and length in lastlengths:  #True if read is a duplicate
		continue
		
	elif lastRead == (chrname, pos) and length not in lastlengths:  #True if read has same start position but different end (length) as previous reads at same start, use this read
		lastlengths.append(length)
		print(line, end="")

	else: #new start position, use this read
		lastRead = (chrname, pos)
		lastlengths=[]
		lastlengths.append(length)
		print(line, end="")