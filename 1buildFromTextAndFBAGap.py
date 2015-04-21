#!/usr/bin/python -tt

import metmodel_current, metmodel_gurobi, wil2metmodelpy, sys
from copy import deepcopy
try:
  reactionfile = sys.argv[1]
except:
  print "Please specify reaction model file (eg: 1buildFromTextAndFBAGap.py my.model.txt)"
name = reactionfile.split('.')[0]
model = metmodel_gurobi.gurobicb()
outprefix = name
reactionfname = reactionfile
exchangefname = '' #''.join((name + '.exchanges.txt'))
biomassfname = ''#''.join((name + 'biomass.txt'))
sourcesfname = ''#''.join((name + 'sources.txt'))
escapesfname = ''#''.join((name + 'escapes.txt'))
epsilon = 0.001

##############################
##  build from text files 
###############################
model.build_from_textfiles(modelfile=reactionfname, exchangesfile=exchangefname, biomassfile=biomassfname, sourcesfile=sourcesfname, escapesfile=escapesfname)
model.convertMet2KeggIDs()

model.addTransports()
## write .wil file
model.writeWil(outprefix + "Working.wil")
model.writeWil(outprefix + "Stable.wil")

################
# build from .wil file
###############
del model
model = metmodel_gurobi.gurobicb()
wilfname = outprefix + "Working.wil"
wil2metmodelpy.wil2metmodel(wilfname)
modelname = wilfname[0:(len(wilfname)-4)]
model.build_from_textfiles(modelfile=modelname + ".reactions", biomassfile=modelname + ".biomass", sourcesfile=modelname+".sources", escapesfile=modelname +".escapes",exchangesfile=exchangefname)

#####################
# add constraints   #
#####################
# reaction id, lower bound, upper bound
#model.set_constraint('R_HEX1',50,50)

model.solve(verbose=False, out=outprefix + "draft")

if model.OBJECTIVE_VALUE < epsilon:
  #######################################################
  # Gap without data base; just add transport reactions #
  #######################################################
  model2 = deepcopy(model)
  model2.fbagapnodb(out=outprefix)
  #############################
  ### Gap fill with data base #
  #############################
  #model3 = deepcopy(model)
  #model3.fbagapdb(out=outprefix)

model.writeECfile(outprefix + "ec1.txt")
metfile = open(outprefix+"metabolites1.txt", "w")
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
