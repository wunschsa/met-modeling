#!/usr/bin/python -tt

import metmodel_current
import metmodel_gurobi
import wil2metmodelpy
from copy import deepcopy

model = metmodel_gurobi.gurobicb()

outprefix = "tfu"
epsilon = 0.001
########################
# build from .wil
########################
wilfname = outprefix + "Working.wil"
wil2metmodelpy.wil2metmodel(wilfname)
modelname = wilfname[0:(len(wilfname)-4)]
model.build_from_textfiles(modelfile=modelname + ".reactions", biomassfile=modelname + ".biomass", sourcesfile=modelname+".sources", escapesfile=modelname +".escapes", exchangesfile="")

#####################
# add constraints   #
#####################
# reaction id, lower bound, upper bound
#model.set_constraint('R_HEX1',50,50)

#######
# FBA #
#######
print "Trying FBA"
model.solve(maps=True, out=outprefix)

model.writeECfile(outprefix + "ec4.txt")
metfile = open(outprefix+"metabolites4.txt", "w")
numspecies = 0
for species in model.SPECIES.keys():
  if not species.endswith("_b"):
    metfile.write("%s\n" % (species))
    numspecies += 1
numreactions = 0
for reactions in model.REACTIONS.keys():
  if (not reactions.startswith("R_SRC")) and (not reactions.startswith("R_ESC")) and (not reactions.startswith("R_EXCH")):
    numreactions += 1
print "Number of metabolites", numspecies
print "Number of reactions", numreactions 
