#!/usr/bin/env python
# coding: utf-8

import math, sys, time, argparse
import mirscanIO as msio

def train(fqueries, bqueries, criteria):
    """
    For each criteria defined in criteria, for each possible output value,
    a score is determined which is the base 2 log of the ratio of the frequency
    of that output value in the foreground vs the frequency of that output
    value in the background. For each criteria, for each output value, the
    foreground and background frequencies are figured such that they include
    pseudocounts that are added according to that feature's 'pseudo' function.

    Output text is in the specified format of a '.matrix' file.
    """

    fscores = criteria.score(fqueries)
    bscores = criteria.score(bqueries)

    # count the occurrence of each value in the foreground and background sets
    fcount = dict()
    bcount = dict()
    for feature in criteria.features:
        fcount[feature] = dict()
        bcount[feature] = dict()
        for v in criteria.features[feature].kl:
            fcount[feature][str(v)] = 0
            bcount[feature][str(v)] = 0
    for cand in fscores:
        for feature in criteria.features:
            fcount[feature][cand[feature]] += 1
    for cand in bscores:
        for feature in criteria.features:
            bcount[feature][cand[feature]] += 1

    output = '# ' + time.asctime(time.localtime()) + '\n'

    # for each criteria, for each value, generate a score and a line for the
    # .matrix output that documents that score and the data that originated it
    for feature in criteria.features:
        output += feature + '\n'

        fcp = criteria.features[feature].pseudo(fcount[feature])
        fcp_total = sum(fcp.values())
        fc_total = sum(fcount[feature].values())

        bcp = criteria.features[feature].pseudo(bcount[feature])
        bcp_total = sum(bcp.values())
        bc_total = sum(bcount[feature].values())

        for v in criteria.features[feature].kl:
            ff = float(fcp[str(v)]) / fcp_total
            bf = float(bcp[str(v)]) / bcp_total
            fraw = round(float(fcount[feature][str(v)]) / fc_total, 3)
            braw = round(float(bcount[feature][str(v)]) / bc_total, 3)
            scr = round(math.log(ff / bf, 2), 3)
            output += '\t'.join(['', str(v), str(scr),\
                                 str(fcount[feature][str(v)]), str(bcount[feature][str(v)]),\
                                 str(fraw), str(braw)]) + '\n'
    return output

parser = argparse.ArgumentParser(description='''MiRscan3 Trainer.
                This script takes a set of foreground (.fam --- the training
                procedure requires mature strand sequences) and background (.fam
                or .fax) miRNA candidates, a mirscan criteria file (.py) with
                rules for their evaluation and outputs a scoring matrix file
                (.matrix) which can be used for scoring additional miRNAs.''',
            epilog='''Paulo Pinto, IEB-WWU, based on:
                http://bartellab.wi.mit.edu/softwareDocs/MiRscan3/Introduction.html''')

parser.add_argument(dest='forefile',
                    help='the foreground training file (.fam)')
parser.add_argument(dest='backfile',
                    help='the background training file (.fam or .fax)')
parser.add_argument(dest='criteriafile',
                    help='the mirscan criteria file (.py)')
parser.add_argument('-o', dest='matrixfile', default='stdout',
                    help='the output scoring matrix file (.matrix) (default: stdout)')

args = parser.parse_args()

# check that args corresponding to filenames have the proper extensions.
if args.forefile.split('.')[-1] != 'fam':
    raise ValueError('Foreground training file must be in \'.fam\' format.')
if args.backfile.split('.')[-1] != 'fam' and args.backfile.split('.')[-1] != 'fax':
    raise ValueError('Background training file must be in \'.fam\' or \'.fax\'  format.')
if args.criteriafile.split('.')[-1] != 'py':
    raise ValueError('Criteria file must be in \'.py\' format.')
if args.matrixfile.lower() != 'stdout' and args.matrixfile.split('.')[-1] != 'matrix':
    raise ValueError('Output matrix file must be in \'.matrix\' format.')

# the foreground set
fqueries = msio.parse_fam(args.forefile)
# the background set
bqueries = msio.parse_query(args.backfile)
# the mirscan criteria
criteria = msio.parse_criteria(args.criteriafile)

with (sys.stdout if args.matrixfile.lower() == 'stdout' else open(args.matrixfile, 'w')) as output:
    output.write(train(fqueries, bqueries, criteria))
