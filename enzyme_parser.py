#!/usr/bin/env python
# enzyme_parser.py
# Written By: Shaun Norris, VCU Bioinformatics MS Candidate
# This will read an ec file downloaded from kegg and extract the necessary information.
# Mostly Interested in the gene information at the time of this writing
import re, sys

filelist = open(sys.argv[1])
filelist = filelist.readlines()

def parse(file):
	keepgoing = True
	infile = open(files)
        outfile = ''.join(( files + "_genes.txt"))
        print outfile
        outfile = open(outfile, "w")
        line = infile.next()
	while keepgoing:
		if re.search('GENES',line):
			print line
			read = 1
			line = re.sub('GENES|\t|\s{2,10}','',line)
			while read > 0:
				if not re.search('[A-Z]{3,4}\:',line) or re.search('DBLINKS',line):
					keepgoing = False
					outfile.close()
					infile.close()
					break
				outfile.write(line + "\n")
				#print line
				line = infile.next()
				line = line.strip()
				
				
		else:
			try:
				line = infile.next()
				line = line.strip()
			except:
				break

for files in filelist:
        files = files.strip()
        print files	
	parse(files)
