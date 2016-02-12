#!/usr/bin/env python2
#####################################################
# Metmodel Pipeline                                 #
# Written By: Shaun Norris VCU Bioinformatics, M.S. #
# Version: 1.0.0  Date: Feb. 2016                   #
#####################################################
""" This is a tool created to generate an in silico model of the metabolic pathways
for a bacterium.  Many of the methods here are moved from, or based off work 
done by Dr. Niti Vanee, Dr. Paul Brooks (VCU), Dr. Steve Fong (VCU) and the whole
effort wouldn't have been possible without the collaborative efforts of 
Stephen Wunsch.
"""
from __future__ import print_function,division
import sys, re, time
import metmodel_current as mc
import metmodel_gurobi as mg
import wil2metmodelpy as w2m
import eq_current
from argparse import ArgumentParser
from copy import deepcopy
from os import path
### Get some args ####
prog = "HOTCHA, the metabolic modeling Pipeline"
ap = ArgumentParser(prog=prog)
ap.add_argument("-i","--in-file",dest="fn",required=True,help="Input the filename/path of the YAML/SEED/Text file you'd like to build from.")
ap.add_argument("-t","--type",dest="ty",required=True,help="Specify the starting filetype. Choices are YAML,SBML,SEED,TXT.")
ap.add_argument("-x","--exchanges",dest="ex",required=True,help="Input the filename contain the exchanges.")
ap.add_argument("-b","--biomass",dest="bmass",required=True,help="Input the filename containg the biomass equation to use.")
ap.add_argument("-d","--experimental-data",dest="exp",default="",required=True,help="Specify experimentally obtained gene/reaction expression data.")
ap.add_argument("-p","--prefix",dest="pfix",default=None,help="Specify a prefix to use for the output files. Default: will be gleaned from the input filename, e.g. 'tfu.init.yaml' becomes prefix of 'tfu'.")
ap.add_argument("-e","--epsilon",dest="ep",default=0.001,type=float,help="Specify an epsilon value to use for the FBA. Default: 0.001")
ap.add_argument("-s","--sources",dest="src",default="",help="Specify a sources file.")
ap.add_argument("-c","--escapes",dest="esc",default="",help="Specify an escapes file.")
ap.add_argument("-r","--biomass-rxn-id",dest="bmrid",default="",help="If using SBML/SEED you MUST specify a biomass reaction ID.")
ap.add_argument("-ndi","--no-data-integration",dest="ndi",action="store_true",default=False,help="Skip integration with experimentally obtained data (proteomic/metabolic/transciptomic/etc.) ")
args = ap.parse_args()
class MetModelPipeline(object):
    def __init__(self,args):
        if args.pfix:
            self.outprefix = args.pfix
        else:
            try:
                self.outprefix = args.fn.split('.')[0]
            except:
                self.outprefix = args.fn
        self.epsilon = args.ep
        self.exchangefname = args.ex
        self.biomassfname = args.bmass
        self.sourcesfname = args.src
        self.escapesfname = args.esc
        self.infile = args.fn
        self.start_type = args.ty
        self.model = mg.gurobicb()
        self.biomassReactionID = args.bmrid
        self.wilfname = self.outprefix + "Working.wil"
        self.modelname = self.wilfname[0:(len(self.wilfname)-4)]

    def build_wil(self):
        if self.start_type.upper() == "YAML":
            self.model.buildfromyaml(self.infile,self.exchangefname)
            self.model.setbiomass(self.biomassfname)
            self.model.convertMet2KeggIDs()
            self.model.addTransports()
            self.model.writeWil(self.outprefix + "Working.wil")
            self.write = True
        elif self.start_type.upper() == "SBML" or self.start_type.upper() == "SEED":
            self.model.read_sbml(self.infile) 
            self.model.convertMet2KeggIDs()  #adds "M_" to metabolite names
            name, reversible, notes, equation = self.model.REACTIONS[biomassReactionID] 
            self.model.add_reaction("R_biomass_target", "BiomassRxn", False, notes, equation) 
            self.model.delete_reaction(biomassReactionID)
            self.model.addTransports()
            self.model.writeWil(self.outprefix + "Working.wil")
            self.write = True
        elif self.start_type.upper() == "TEXT" or self.start_type.upper() == "TXT":
            self.model.build_from_textfiles(modelfile=self.infile, exchangesfile=self.exchangefname, biomassfile=self.biomassfname, sourcesfile=self.sourcesfname, escapesfile=self.escapesfname)
            self.model.convertMet2KeggIDs()
            self.model.addTransports()
            self.model.writeWil(self.outprefix + "Working.wil")
            self.model.writeWil(self.outprefix + "Stable.wil")
            self.write = False
        del self.model
        self.model = mg.gurobicb()
        w2m.wil2metmodel(self.wilfname)
        self.model.build_from_textfiles(modelfile=self.modelname + ".reactions", biomassfile=self.modelname + ".biomass", sourcesfile=self.modelname+".sources", escapesfile=self.modelname +".escapes",exchangesfile="")
        if self.write:
            self.model.writeWil(self.outprefix + "Working.wil")
            self.model.writeWil(self.outprefix + "Stable.wil")

    def solve(self,verbose = False,withdb = False):
        self.model.solve(verbose=verbose,out=self.outprefix + 'draft')

        if self.model.OBJECTIVE_VALUE < self.epsilon:
            model2 = deepcopy(self.model)
        if withdb:
            model2.fbagapdb(out=self.outprefix)
        else:
            model2.fbagapnodb(out=self.outprefix)

        self.model.writeECfile(self.outprefix  + "ec1.txt")
        metfile = open(self.outprefix + "metabolites1.txt","w")
        numspecies = 0
        for species in self.model.SPECIES.keys():
            if not species.endswith("_b"):
                metfile.write("%s\n" % (species))
                numspecies += 1
        numreactions = 0
        for reactions in self.model.REACTIONS.keys():
            if (not reactions.startswith("R_SRC")) and (not reactions.startswith("R_ESC")) and (not reactions.startswith("R_EXCH")):
                numreactions += 1
        del model2
        del self.model
        return numspecies,numreactions

    def gapfill(self):
        self.model = mg.gurobicb()
        w2m.wil2metmodel(self.wilfname)
        self.model.build_from_textfiles(modelfile=self.modelname + ".reactions", biomassfile=self.modelname + ".biomass", sourcesfile=self.modelname+".sources", escapesfile=self.modelname +".escapes", exchangesfile="")
        try:
            print ("Reading reactions to add")
            gapfile = open(self.outprefix + ".gap.xls", "r")
            reactionsToAdd = []
            source = 0
            escape = 0
            reaction = 0
            for line in gapfile.readlines():
                if "added sources:" in line:
                    source = 1
                    continue
                if "added escapes:" in line:
                    escape = 1
                    source = 0
                    continue
                if "added reactions:" in line:
                    reaction = 1
                    escape = 0
                    continue
                if reaction == 1:
                    reactionsToAdd.append(line.rstrip())
                if escape == 1:
                    escapeToAdd = line.rstrip()
                    self.model.ESCAPES.append(escapeToAdd)
                if source == 1:
                    sourceToAdd = line.rstrip()
                    self.model.SOURCES.append(sourceToAdd)
                self.model.addReactionsFromDB(reactionsToAdd)  
                self.model.set_sources(model.SOURCES)
                self.model.set_escapes(model.ESCAPES)
        except:
            print ("No .gap.xls file")
        self.model.writeWil(filename=self.outprefix+"Working.wil")

    def fba(self):
        print ("Trying FBA")
        if self.model.OBJECTIVE_VALUE < self.epsilon:
            print ("Flux was zero, gap filling...")
            model3 = deepcopy(model)
            model3.fbagapdb(out=self.outprefix)
        else:
            print ("Flux was positive.")
            gapfile = open(self.outprefix + ".gap.xls", "w")
            gapfile.write("added sources:\nadded escapes:\nadded reactions:\n")
            gapfile.close()
        try:
            del model3
        except:
            pass
        self.model.writeECfile(self.outprefix + "ec2.txt")
        metfile = open(self.outprefix + "metabolites2.txt","w")
        numspecies = 0
        for species in self.model.SPECIES.keys():
            if not species.endswith("_b"):
                metfile.write("%s\n" % (species))
                numspecies += 1
        numreactions = 0
        for reactions in self.model.REACTIONS.keys():
            if (not reactions.startswith("R_SRC")) and (not reactions.startswith("R_ESC")) and (not reactions.startswith("R_EXCH")):
                numreactions += 1
        del self.model
        return numspecies,numreactions

    def data_integration(self,gammas = (-0.1,0.1),exprfile = args.exp,forceBiomass = True, verbose = False):
        """
        Ruppin  analysis of syw using the datasets from Tetu et al., 2009 about phosphate stress.
        """
        gprfname = self.infile + '.gpr'
        def makeCallsDict( allRxns, hirxns, lorxns ):
            calls = {}
            for r in hirxns:
                assert r in allRxns and r not in lorxns
            for r in lorxns:
                assert r in allRxns and r not in hirxns
            for r in allRxns:
                if r in hirxns:
                    calls[r] = "On"
                elif r in lorxns:
                    calls[r] = "Off"
                else: 
                    calls[r] = "No call"
            for r in allRxns:
                assert not calls[r] == None
            return calls
        self.model = mg.gurobicb()
        w2m.wil2metmodel(self.wilfname)
        self.model.build_from_textfiles(modelfile=self.modelname + ".reactions", biomassfile=self.modelname + ".biomass", sourcesfile=self.modelname+".sources", escapesfile=self.modelname +".escapes", exchangesfile="")
        if verbose:
            rq = False
        else:
            rq = True
        self.model.gpr2(gprfname,readquiet=rq)

        fluxes = {}
        calls = {}

        self.model.solve(verbose=verbose)
        fluxes['fba'] = deepcopy(self.model.REACTION2FLUXVALUE)

        self.model.ruppin(exprfile,gammas,forceBiomass=forceBiomass)

        calls[self.outprefix] = makeCallsDict(self.model.REACTIONS,self.model.hirxnset,self.model.lowrxnset)
        fluxes[self.outprefix] = deepcopy(self.model.REACTION2FLUXVALUE)
        cache = {}
        for reaction in self.model.REACTIONS:
            name, reversible, notes, equation = self.model.REACTIONS[reaction]
            reactionequation = eq_current.makestring(equation, reversible)
            confidence,gpr = '?','?'
            holder = {'pathways':{}, 'ecs':{}}
            ref, prr = '.', '.'
            for note in notes:
                if 'SUBSYSTEM: ' in note:
                    holder['pathways'][note[11:]] = 1
                if 'EC: ' in note:
                    holder['ecs'][note[4:]] = 1
                if 'CONFIDENCE: ' in note:
                    confidence = note[12:]
                if 'GPR: ' in note:
                    gpr = note[5:]
                if 'Protein_reaction_relation: ' in note:
                    prr = note[note.find(' == ') + 4:]
                if 'PMID: ' in note: # and not 'review' in note and not 'related_organism' in note 
                    ref = ref + note.split(',')[0][6:] + ' '
            if not ref == '.':
                ref = 'PMIDs: ' + ref[1:-1]

            grr = prr
            for i in prr.split():
                if '(' in i:
                    i = i[1:]
                if ')' in i:
                    i = i[:-1]
                if i in self.model.PROTEIN2GENE:
                    grr = grr.replace(i, model.PROTEIN2GENE[i])
                    
            #add to cache for printing...                                                                       
            for pathwayname in holder['pathways']:
                for ec in holder['ecs']:
                    if not pathwayname in cache:
                        cache[pathwayname] = {}
                    if not ec in cache[pathwayname]:
                        cache[pathwayname][ec] = {}

                    cache[pathwayname][ec][ ('\t').join((reaction, name, str(reversible), pathwayname, ec, fluxes['fba'][reaction], fluxes[self.outprefix][reaction], reactionequation)) ] = 1
        paths = cache.keys()
        paths.sort()
        outfile = open(self.outprefix+"integrated.xls","w")

        for path in paths:
            ecs = cache[path].keys()
            ecs.sort()
            for ec in ecs:
                for r in cache[path][ec]:
                    outfile.write("%s\n" % r)
        del self.model

    def mapFlux(self):
        self.model = mg.gurobicb()
        
        w2m.wil2metmodel(self.wilfname)
        self.model.build_from_textfiles(modelfile=self.modelname + ".reactions", biomassfile=self.modelname + ".biomass", sourcesfile=self.modelname+".sources", escapesfile=self.modelname +".escapes",exchangesfile="")
        self.model.solve(maps=True,out=self.outprefix)
        self.model.writeECfile(self.outprefix + "ec4.txt")
        metfile = open(self.outprefix+"metabolites4.txt", "w")
        numspecies = 0
        for species in self.model.SPECIES.keys():
          if not species.endswith("_b"):
            metfile.write("%s\n" % (species))
            numspecies += 1
        numreactions = 0
        for reactions in self.model.REACTIONS.keys():
          if (not reactions.startswith("R_SRC")) and (not reactions.startswith("R_ESC")) and (not reactions.startswith("R_EXCH")):
            numreactions += 1
        return numspecies,numreactions


def main(args):
    model = MetModelPipeline(args)
    model.build_wil()
    num_metabolites,num_reactions = model.solve()
    print ("Step 1. complete with, %i reactions and %i metabolies." % (num_metabolites,num_reactions))
    model.gapfill()
    num_metabolites,num_reactions = model.fba()
    print ("Step 2. Gap Fill & FBA complete with, %i reactions and %i metabolites."% (num_metabolites,num_reactions))
    if not args.ndi:
        model.data_integration(exprfile=args.exp)
        print ("Step 3. Data Integration is now complete.")
    num_metabolites,num_reactions = model.mapFlux()
    print ("Final step is now complete, model is ready.")
    print ("Final model contains %i reactions and %i metabolites" % (num_metabolites,num_reactions))

if __name__ == '__main__':
    if len(sys.argv) < 8:
        print ("Insufficient arguments provided...")
        ap.print_help()
        exit(1)
    else:
        if (path.isfile(args.fn) and path.isfile(args.ex) and path.isfile(args.bmass)):
            main(args)
        else:
            files = [args.fn,args.ex,args.bmass]
            badfiles = [path.isfile(args.fn), path.isfile(args.ex), path.isfile(args.bmass)]
            if args.src:
                files.append(args.src)
                badfiles.append(path.isfile(args.src))
            if args.esc:
                files.append(args.esc)
                badfiles.append(path.isfile(args.esc))
            for each in range(0,len(files)):
                try:
                    print ("Could not find: %s" % files.pop(badfiles.index(False)))
                    badfiles.pop(badfiles.index(False))
                except:
                    pass
            print ("Please check filename(s) & paths then try again.")
            exit(1)
