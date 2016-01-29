#!/usr/bin/env python2
# Asgard Pipeline #
# Written by: Shaun Norris
# /usr/global/blp/bin/asgard -i tfu.fasta -p blastn -n 20 -d /gpfs_fs/data/refdb/asgardDB/UniRef100 -d /gpfs_fs/data/refdb/asgardDB/KEGG -f /gpfs_fs/data/refdb/asgardDB/uniref100.fasta.gz -f /gpfs_fs/data/refdb/asgardDB/genes.pep.gz -l /gpfs_fs/data/refdb/asgardDB
from __future__ import print_function,division
import sys,re
import subprocess as sp

class Asgard(object):
    def __init__(self,fasta,blast = "blastn",nodes = 10,threads = 4,dblist=["/gpfs_fs/data/refdb/asgardDB/UniRef100","/gpfs_fs/data/refdb/asgardDB/KEGG"],fa_dblist=["/gpfs_fs/data/refdb/asgardDB/uniref100.fasta.gz","/gpfs_fs/data/refdb/asgardDB/genes.pep.gz"],mapping_fp="/gpfs_fs/data/refdb/asgardDB"):
        self.infile = fasta
        self.blast = blast
        self.nodes = nodes
        self.threads = 4
        self.dblist = "-d ".join(dblist)
        self.fa_dblist = "-f ".join(fa_dblist)
        self.mapping_fp = mapping_fp

    def run(self):
        sp.Popen(['asgard','-i',self.infile,'-p',self.blast,'-n',self.nodes,self.db_list,self.fa_dblist,'-l',mapping_fp]).wait()

tfu = Asgard('tfu.fasta')
tfu.run()