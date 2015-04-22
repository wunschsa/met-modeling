#!/usr/bin/env python
# enzyme_parser.py
# Written By: Shaun Norris, VCU Bioinformatics MS Candidate
# This will read an ec file downloaded from kegg and extract the necessary information.
# Mostly Interested in the gene information at the time of this writing
import re, sys,psycopg2
from collections import defaultdict

filelist = open(sys.argv[1])
filelist = filelist.readlines()

def parse(file):
	defdict = defaultdict(defaultdict)
	keepgoing = True
	#outfile = ''.join(( files + "_genes.txt"))
	outfile = open("databaseready.tsv", "a")
	enzyme = re.sub('ec:|\;','',re.split('\/',file)[-1])
	#outfile.write(re.sub('ec:','',re.split('\/',file)[-1]))
	infile = open(files)
	line = infile.next()
	while keepgoing:
		if re.match('NAME',line):
			defdict[enzyme]['name'] = re.sub('NAME\s*|\;','',line)
			#outfile.write('\t%s' % re.sub('NAME\s*|\;','',line))
		if re.match('REACTION',line):
			defdict[enzyme]['rxn'] = re.sub('REACTION\s*|\;','',line)
			#outfile.write('\t%s' % re.sub('REACTION\s*','',line))
		if re.match('ALL_REAC',line):
			defdict[enzyme]['rxnids'] = re.sub('ALL_REAC\s*|\;','',line)
			#outfile.write('\t%s' % re.sub('ALL_REAC\s*','',line))
		if re.search('GENES',line):
			read = 1
			defdict[enzyme]['genes'] = re.sub('GENES|\t|\s{2,10}|\;','',line)
			line = re.sub('GENES|\t|\s{2,10}','',line)
			#outfile.write('\t' + line + " ")
			try:
				line = infile.next()
				line = line.strip()
			except:
				break
			while read > 0:
				if not re.search('.*[A-Z]{3,4}\:',line) or re.search('DBLINKS',line):
					keepgoing = False
					#outfile.write("\n")
					#outfile.close()
					infile.close()
					break
				defdict[enzyme]['genes'] += " " + re.sub('GENES|\t|\s{2,10}|\;','',line)
				#outfile.write(line + " ")
				line = infile.next()
				line = line.strip()
				
				
		else:
			try:
				line = infile.next()
				line = line.strip()
			except:
				#outfile.write("\n")
				break
	for key in defdict:
		if 'name' not in defdict[key]:
			defdict[key]['name'] = "None"
		if 'rxn' not in defdict[key]:
			defdict[key]['rxn'] = "None"
		if 'rxnids' not in defdict[key]:
			defdict[key]['rxnids'] = "None"
		if 'genes' not in defdict[key]:
			defdict[key]['genes'] = "None"
		#outfile.write('%s\t%s\t%s\t%s\t%s\n' % key,defdict[key]['name'],defdict[key]['rxn'],defdict[key]['rxnids'],defdict[key]['genes'])
		#except:
		outfile.write("%s\t%s\t%s\t%s\t%s\n" % (key,defdict[key]['name'],defdict[key]['rxn'],defdict[key]['rxnids'],defdict[key]['genes']))
		#outfile.write('%s\t%s\t%s\t%s\t%s\n' % key,defdict[key]['name'],defdict[key]['rxn'],defdict[key]['rxnids'],"None")
	outfile.close()

for files in filelist:
	files = files.strip()
	parse(files)
try:
    conn = psycopg2.connect("dbname='reactions' user='snorris' host='localhost'")#" password='testing!@3'") # open that database connection  #other user snake password pyth0n! dbname='reactions'
except:
    print "Unable to connect to database"
    quit()
cur = conn.cursor()
cur.copy_from(open("databaseready.tsv"),'kegg_enzymes')
conn.commit()
conn.close()
