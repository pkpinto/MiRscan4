#!/usr/bin/env python
# coding: utf-8

# This script applies a set of scoring criteria to a set of candidate miRNA hairpins
# and outputs the resulting score sheet.

import sys, argparse
import mirscanModule



parser = argparse.ArgumentParser(description='MiRscan3 Execute',
            epilog='Paulo Pinto, IEB-WWU, based on:\nhttp://bartellab.wi.mit.edu/softwareDocs/MiRscan3/Introduction.html')

parser.add_argument(dest='queryFile',
                    help='the query file (.train or .fax)')
parser.add_argument(dest='criteriaFile',
                    help='the mirscan criteria file (.py)')
parser.add_argument(dest='matrixFile',
                    help='the scoring matrix file (.matrix)')

parser.add_argument('-o', dest='scoreFile', default='stdout',
                    help='the output score sheet file (.scr) (default: stdout)')

args = parser.parse_args()

# check that args corresponding to filenames have the proper extensions.
if args.queryFile.split('.')[-1]!='train' and args.queryFile.split('.')[-1]!='fax':
    raise ValueError('query file must be ".train" or ".fax" format')
if args.criteriaFile.split('.')[-1]!='py':
    raise ValueError('criteria file must be formatted for python (".py")')
if args.matrixFile.split('.')[-1]!='matrix':
    raise ValueError('scoring matrix must be ".matrix" format')
if args.scoreFile.lower()!='stdout' and args.scoreFile.split('.')[-1]!='scr':
    raise ValueError('outfile must be ".scr".')


# get the mirscan criteria dictionary, 'fdict'
f = open(args.criteriaFile)
mirscanText = f.read()
f.close()
mirscanDict = dict()
exec mirscanText in mirscanDict
fdict = mirscanDict['fdict']


# get the scoring matrix and queries
ms = mirscanModule.get_ms(args.matrixFile)
queryList = mirscanModule.get_queries(args.queryFile)


with (sys.stdout if args.scoreFile.lower()=='stdout' else open(args.scoreFile,'w')) as output:

    # score each candidate and print out the results in ".scr" specified format.
    for candScore in mirscanDict['mirscan'](queryList,ms):
        data = []
        data.append(candScore['name'])
        data.append('totscore '+str(candScore['totscore']))
        # print out loc values for each species
        locKeys = filter(lambda k: k[:4]=='loc_', candScore.keys())
        for k in locKeys: data.append(k+' '+str(candScore[k]))
        for k in fdict.keys(): data.append(k+' '+str(candScore[k]))
        output.write(' '.join(data)+'\n')
