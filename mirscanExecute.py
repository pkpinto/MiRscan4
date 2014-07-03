# This script applies a set of scoring criteria to a set of candidate miRNA hairpins
# and outputs the resulting score sheet.
#
# argv[1] is the query file (.train or .fax)
# argv[2] is the mirscan criteria file (.py)
# argv[3] is the scoring matrix file (.matrix)
# argv[4] is the output file (.scr)

import sys, mirscanModule

queryFile = sys.argv[1]
mirscanFile = sys.argv[2]
matrixFile = sys.argv[3]
outfile = sys.argv[4]

if queryFile.split('.')[-1]!='train' and queryFile.split('.')[-1]!='fax':
    raise ValueError('query file must be ".train" or ".fax" format')
if mirscanFile.split('.')[-1]!='py':
    raise ValueError('mirscan criteria file must be formatted for python (".py")')
if matrixFile.split('.')[-1]!='matrix':
    raise ValueError('scoring matrix must be ".matrix" format')
if outfile.lower()!='stdout' and outfile[-4:]!='.scr':
    raise ValueError("outfile name should end with '.scr'.")


# get the mirscan criteria dictionary, 'fdict'
f = open(mirscanFile)
mirscanText = f.read()
f.close()
mirscanDict = dict()
exec mirscanText in mirscanDict
fdict = mirscanDict['fdict']


# get the scoring matrix and queries
ms = mirscanModule.get_ms(matrixFile)
queryList = mirscanModule.get_queries(queryFile)



if outfile.lower()=='stdout': out_to = sys.stdout
else: out_to = open(outfile,'w')



# score each candidate and print out the results in ".scr" specified format.
for candScore in mirscanDict['mirscan'](queryList,ms):
    data = []
    data.append(candScore['name'])
    data.append('totscore '+str(candScore['totscore']))
    # print out loc values for each species
    locKeys = filter(lambda k: k[:4]=='loc_', candScore.keys())
    for k in locKeys: data.append(k+' '+str(candScore[k]))
    for k in fdict.keys(): data.append(k+' '+str(candScore[k]))
    out_to.write(' '.join(data)+'\n')

out_to.close()
