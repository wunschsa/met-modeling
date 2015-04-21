import re, psycopg2
try:
    conn = psycopg2.connect("dbname='reactions' user='creator' host='localhost' password='testing!@#'") # open that database connection  #other user snake password pyth0n! dbname='reactions'
except:
    print "Unable to connect to database"
    quit()
cur = conn.cursor() #set up the psql cursor

modelfile = open('/home/snorris/Ubuntu One/Thesis Project/Models/Pseudomonas_aeruginosa/smallPaeruginosa.csv','r')
model = modelfile.readlines()
modelfile.close()
outfile = open('/home/snorris/Ubuntu One/Thesis Project/Models/Pseudomonas_aeruginosa/Paeruginosa.wil','w')
for line in model:
    if re.match('[0-9]',line):
        line = line.split('\t')
        name = line[2]
        EC = line[9]
        rxns = re.split(' ',line[3])
        newrxn = ''
        kegg_rxnID = ''
        true = 'True'
        null = 'null'
        for each_rxn in range(0,len(rxns)):
            #print rxns[each_rxn]
            if re.search('\(|\)|\+|\>|\<|\-\-|\=|\:|\s|\[|\]',rxns[each_rxn]):
                newrxn += rxns[each_rxn]
           #     print newrxn
            else:
                cur.execute(("SELECT keggcpd from kegg_compounds where to_tsvector(names) @@ to_tsquery('%s') ORDER BY keggcpd;") % rxns[each_rxn])
                newrxn += re.sub('\'|\(|\)|\,','',str(cur.fetchone()))
        outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (name, EC, kegg_rxnID, true, name, newrxn, null, null, null))
            
    