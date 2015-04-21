#!/usr/bin/python -tt

import metmodel_current, metmodel_gurobi
import wil2metmodelpy
from copy import deepcopy
import sys
try:
  name = str(sys.argv[1]).split('.')[0]
except:
  print "Please input the name or abbreviation of the organism"
  sys.exit(2)
model = metmodel_gurobi.gurobicb()

outprefix = name
epsilon = 0.001
########################
# build from .wil
########################
wilfname = outprefix + "Working.wil"
wil2metmodelpy.wil2metmodel(wilfname)
modelname = wilfname[0:(len(wilfname)-4)]
model.build_from_textfiles(modelfile=modelname + ".reactions", biomassfile=modelname + ".biomass", sourcesfile=modelname+".sources", escapesfile=modelname +".escapes", exchangesfile="")

# read in curated transports and reactions to add
try:
  print "Reading reactions to add"
  gapfile = open(outprefix + ".gap.xls", "r")
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
      model.ESCAPES.append(escapeToAdd)
    if source == 1:
      sourceToAdd = line.rstrip()
      model.SOURCES.append(sourceToAdd)
  model.addReactionsFromDB(reactionsToAdd)  
  model.set_sources(model.SOURCES)
  model.set_escapes(model.ESCAPES)
except:
  print "No .gap.xls file"



# write new .wil file
model.writeWil(filename=outprefix+"Working.wil")

#####################
# add constraints   #
#####################
# reaction id, lower bound, upper bound
#model.set_constraint('R_HEX1',50,50)

#######
# FBA #
#######
print "Trying FBA"
model.solve(verbose=False, out=outprefix)

# if still no flux, do gap analysis; otherwise, write new model to .wil file
if model.OBJECTIVE_VALUE < epsilon:
  print "Flux was zero, gap filling"
  ########################################################
  ## Gap without data base; just add transport reactions #
  ########################################################
  #model2 = deepcopy(model)
  #model2.fbagapnodb(out=outprefix)
  ############################
  ## Gap fill with data base #
  ############################
  model3 = deepcopy(model)
  model3.fbagapdb(out=outprefix)
else:
  print "Flux was positive"
  gapfile = open(outprefix + ".gap.xls", "w")
  gapfile.write("added sources:\nadded escapes:\nadded reactions:\n")
  gapfile.close()
    
model.writeECfile(outprefix + "ec2.txt")
metfile = open(outprefix+"metabolites2.txt", "w")
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
