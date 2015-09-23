#!/usr/bin/python -tt

import metmodel_current
import metmodel_gurobi
import wil2metmodelpy
from copy import deepcopy

model = metmodel_gurobi.gurobicb()
outprefix = "Efaulk"
exchangefname = ""
biomassReactionID = "bio23103"
sbmlFname = "models/Efaulk_RAST.xml"
epsilon = 0.001

#################################
# Read from SBML
#################################

model.read_sbml(sbmlFname) 
model.convertMet2KeggIDs()  #adds "M_" to metabolite names
name, reversible, notes, equation = model.REACTIONS[biomassReactionID] 
model.add_reaction("R_biomass_target", "BiomassRxn", False, notes, equation) 
model.delete_reaction(biomassReactionID)

model.addTransports()
model.writeWil(outprefix + "Working.wil")

################
## build from .wil file
################         
del model             
model = metmodel_gurobi.gurobicb()
wilfname = outprefix + "Working.wil"
wil2metmodelpy.wil2metmodel(wilfname)
modelname = wilfname[0:(len(wilfname)-4)]
model.build_from_textfiles(modelfile=modelname + ".reactions", biomassfile=modelname + ".biomass", sourcesfile=modelname+".sources", escapesfile=modelname +".escapes",exchangesfile=exchangefname)

# write .wil file with exchanges
model.writeWil(outprefix + "Working.wil")
model.writeWil(outprefix + "Stable.wil")

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
