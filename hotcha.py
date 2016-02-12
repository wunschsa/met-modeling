#!/usr/bin/env python2
# Theoretically this should be compatible with python 3 as well. #
#############################################################
# Hotcha - The Asgard Pipeline                              #
# Written by: Shaun Norris, VCU M.S. Bioinformatics         #
# Contact: jpbrooks [at] vcu.edu OR norrissw [at] vcu.edu   #
# Current version 0.0.1a - Interal use only                 #
#############################################################
from __future__ import print_function,division
import sys,re
import subprocess as sp

class Asgard(object):
    def __init__(self,fasta,blast = "blastn",nodes = 10,threads = 4,dblist=["/gpfs_fs/data/refdb/asgardDB/UniRef100","/gpfs_fs/data/refdb/asgardDB/KEGG"],fa_dblist=["/gpfs_fs/data/refdb/asgardDB/uniref100.fasta.gz","/gpfs_fs/data/refdb/asgardDB/genes.pep.gz"],mapping_fp="/gpfs_fs/data/refdb/asgardDB"):
        self.infile = fasta
        self.blast = blast
        self.nodes = nodes
        self.threads = 4
        self.dblist = [] 
        self.fa_dblist = []
        self.mapping_fp = mapping_fp
        for db_path in dblist:
            self.dblist.append("-d")
            self.dblist.append(db_path)
        for fa_dbpath in fa_dblist:
            self.fa_dblist.append("-f")
            self.dblist.append(fa_dblist)


    def run(self):
        """ 
        This will actually run the pipeline
        Example of the general usage of asgard:
        asgard -i tfu.fasta -p blastn -n 20 -d UniRef100 -d KEGG -f uniref100.fasta.gz -f genes.pep.gz -l asgardDB
        """
        sp.Popen(['asgard','-i',self.infile,'-p',self.blast,'-n',self.nodes,self.dblist,self.fa_dblist,'-l',self.mapping_fp]).wait()

if __name__ == '__main__':
    tfu = Asgard('tfu.fasta') # example usage of creating an instance of the Asgard Object
    tfu.run() # run the pipeline