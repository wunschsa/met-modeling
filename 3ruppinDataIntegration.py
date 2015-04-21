#!/usr/bin/python -tt

# Ruppin  analysis of syw using the datasets from Tetu et al., 2009 about phosphate stress.
import sys
try:
	name = str(sys.argv[1]).split('.')[0]
except:
	print "Please input the name or abbreviation of the organism"
	sys.exit(2)
outprefix = name
gprfname = ''.join((name + '.gpr'))

gammas = (-0.1, 0.1)	# correlates to type I error of alpha = 0.00005

#defaultConstraints = {'R_EXCH_glc_DASH_D_e': (-100, 0) }

experimentDict = { 
	'TFU_Proteomics_Cb2D' : { 'exprfile': 'models/TFU_Proteomics_Cb2D.txt' }
	# Commented following files on 05_09_2012 by Niti Vanee
#	'F41_M9_7DS521' : { 'exprfile':'../../../Datasets/RDE2_X3339_FlightvsGround_F41_M9_7DS521.txt' },
#	'F423_M9_7DS525' : { 'exprfile':'../../../Datasets/RDE2_X3339_FlightvsGround_F423_M9_7DS525.txt'},
#	'F54_M9_7DS522' : { 'exprfile':'../../../Datasets/RDE2_X3339_FlightvsGround_F54_M9_7DS522.txt' },
#	'G412_M9_7DS523_only2subs' : { 'exprfile':'../../../Datasets/RDE2_X3339_FlightvsGround_G412_M9_7DS523_only2subs.txt' },
#	'G43_M9_7DS524' : { 'exprfile':'../../../Datasets/RDE2_X3339_FlightvsGround_G43_M9_7DS524.txt'},
#	'G514_M9_7DS526' : { 'exprfile':'../../../Datasets/RDE2_X3339_FlightvsGround_G514_M9_7DS526.txt' }
		}
	
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
	
# build model:
import metmodel_current
import wil2metmodelpy
import metmodel_gurobi
import eq_current

from copy import deepcopy
model = metmodel_gurobi.gurobicb()

wilfname = outprefix + "Working.wil"
wil2metmodelpy.wil2metmodel(wilfname)
modelname = wilfname[0:(len(wilfname)-4)]
model.build_from_textfiles(modelfile=modelname + ".reactions", biomassfile=modelname + ".biomass", sourcesfile=modelname+".sources", escapesfile=modelname +".escapes", exchangesfile="")
model.gpr2(gprfname, readquiet=True)

	
fluxes = {}
calls = {}

#Run fba experiment:
model.solve(verbose=False)
fluxes['fba'] = deepcopy(model.REACTION2FLUXVALUE)
print fluxes['fba']['R_biomass_target']

#Run ruppin experiments:
for x in experimentDict:
	
	print x
	#set constraints:
	#constraints = experimentDict[x].get('constraints', {})
	#for c in constraints:
	#	model.set_constraint(c[0], c[1], c[2])
	
	exprfile = experimentDict[x]['exprfile']
	model.ruppin(exprfile, gammas, forceBiomass=True)
	calls[x] = makeCallsDict( model.REACTIONS, model.hirxnset, model.lowrxnset )
	fluxes[x] = deepcopy(model.REACTION2FLUXVALUE)
	
	#reset default constraints:
	#for c in constraints:
	#	(lb, ub) = defaultConstraints[c[0]]
	#	model.set_constraint(c[0], lb, ub)
		
#model.print_flux_xls( fluxes, 'Ruppins_SpaceFlight_05_09_2012.xls', calls )
	cache = {}
	for reaction in model.REACTIONS:
		name, reversible, notes, equation = model.REACTIONS[reaction]
		reactionequation = eq_current.makestring(equation, reversible)
		#search for any ec numbers and pathways in reaction notes; it IS possible for there to be > 1 ec or pathway for a given reaction
		confidence, gpr = '?', '?'
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
		
		#construct grr
		grr = prr
		for i in prr.split():
			if '(' in i:
				i = i[1:]
			if ')' in i:
				i = i[:-1]
			if i in model.PROTEIN2GENE:
				grr = grr.replace(i, model.PROTEIN2GENE[i])
				
		#add to cache for printing...																		
		for pathwayname in holder['pathways']:
			for ec in holder['ecs']:
				if not pathwayname in cache:
					cache[pathwayname] = {}
				if not ec in cache[pathwayname]:
					cache[pathwayname][ec] = {}
					
				#if this is printing the results of 'solve', then get flux activities...
				#reactionID	name	rev	pathway	ec	equation	confidence	gpr
				#if showfluxvalues and self.REACTION2FLUXVALUE.get(reaction, '.') != '0':
				cache[pathwayname][ec][ ('\t').join((reaction, name, str(reversible), pathwayname, ec, fluxes['fba'][reaction], fluxes[x][reaction], reactionequation)) ] = 1
	
	paths = cache.keys()
	paths.sort()		
	outfi = open(outprefix+"Ruppin.xls", 'w')
		
	for path in paths:
		ecs = cache[path].keys()
		ecs.sort()
		for ec in ecs:
			for r in cache[path][ec]:
				print >>outfi, r
		#print >>outfi, '\n'
			
