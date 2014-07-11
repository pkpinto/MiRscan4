#!/usr/bin/env python
# coding: utf-8

import sys, argparse
import mirscanIO as msio

parser = argparse.ArgumentParser(description='''MiRscan3 Scorer.
                Given a query file listing miRNA candidates (.fam or .fax),
                this script scores them using an (also given) mirscan
                criteria file (.py) and scoring matrix file (.matrix). The
                output is done to a score sheet file (.scr).''',
            epilog='''Paulo Pinto, IEB-WWU, based on:
                http://bartellab.wi.mit.edu/softwareDocs/MiRscan3/Introduction.html''')

parser.add_argument(dest='queryfile',
                    help='the query file (.fam or .fax)')
parser.add_argument(dest='criteriafile',
                    help='the mirscan criteria file (.py)')
parser.add_argument(dest='matrixfile',
                    help='the scoring matrix file (.matrix)')
parser.add_argument('-o', dest='scorefile', default='stdout',
                    help='the output score sheet file (.scr) (default: stdout)')

args = parser.parse_args()

# check that args corresponding to filenames have the proper extensions.
if args.queryfile.split('.')[-1] != 'fam' and args.queryfile.split('.')[-1] != 'fax':
    raise ValueError('Query file must be in \'.fam\' or \'.fax\' format.')
if args.criteriafile.split('.')[-1] != 'py':
    raise ValueError('Criteria file must be in \'.py\' format.')
if args.matrixfile.split('.')[-1] != 'matrix':
    raise ValueError('Matrix file must be in \'.matrix\' format.')
if args.scorefile.lower() != 'stdout' and args.scorefile.split('.')[-1] != 'scr':
    raise ValueError('Output score sheet file must be in \'.scr\' format.')


candidates = msio.parse_query(args.queryfile)
mse = msio.parse_criteria(args.criteriafile)
matrix = msio.parse_matrix(args.matrixfile)

with (sys.stdout if args.scorefile.lower() == 'stdout' else open(args.scorefile, 'w')) as output:

    # score each candidate and print out the results in ".scr" specified format.
    for candScore in mse['mirscan'](candidates, matrix):
        data = []
        data.append(candScore['name'])
        data.append('totscore '+str(candScore['totscore']))
        # print out loc values for each species
        locKeys = filter(lambda k: k[:4]=='loc_', candScore.keys())
        for k in locKeys: data.append(k+' '+str(candScore[k]))
        for k in mse['fdict'].keys(): data.append(k+' '+str(candScore[k]))
        output.write(' '.join(data)+'\n')
