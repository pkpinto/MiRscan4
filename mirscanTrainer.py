#!/usr/bin/env python
# coding: utf-8

import math, sys, time, argparse
import mirscanIO as msio

def train(fqueries, bqueries, criteria):
    """
    Given a set of foreground (fqueries, must include information where the
    mature strands are located in the hairpin) and background candidates
    (bqueries), this method first computes the raw scores for all features
    described in the critera. Then, separately for the foreground and the
    background set, the frequency of all possible values for each feature is
    computed.

    For the ouput, each feature is considered in turn. The distribution over
    possible values for a particular feature is smoothed by calling the
    criteria.features[feature].smooth function. The final score for each value
    is given by the logarithm base 2 of the ratio between the (smoothed)
    frequency measured for the foreground set to the (smoothed) frequency in the
    background set.

    Further output includes, per feature value, of the counts made in the
    foreground and the background sets and the corresponding unsmoothed
    frequencies.

    The function returns a string formated in the .matrix format.
    """
    fscores = criteria.score(fqueries)
    bscores = criteria.score(bqueries)

    # count the occurrence of each value in the foreground and background sets
    fcount = dict()
    bcount = dict()
    for feature in criteria.features:
        fcount[feature] = {b: 0 for b in criteria.features[feature].bins}
        bcount[feature] = {b: 0 for b in criteria.features[feature].bins}
    for cand_score in fscores:
        for feature in criteria.features:
            fcount[feature][criteria.features[feature].pick_bin(cand_score[feature])] += 1
    for cand_score in bscores:
        for feature in criteria.features:
            bcount[feature][criteria.features[feature].pick_bin(cand_score[feature])] += 1

    output = '# ' + time.asctime(time.localtime()) + '\n'

    for feature in criteria.features:
        output += feature + '\n'

        # sum counts for feature values
        ffeature_sum = sum(fcount[feature].values())
        bfeature_sum = sum(bcount[feature].values())
        # smooth the count distribution
        ffeature_smooth = criteria.features[feature].smooth(fcount[feature])
        bfeature_smooth = criteria.features[feature].smooth(bcount[feature])
        # and sum the values
        ffeature_smooth_sum = sum(ffeature_smooth.values())
        bfeature_smooth_sum = sum(bfeature_smooth.values())

        for b in criteria.features[feature].bins:
            ffreq_smooth = 1.0 * ffeature_smooth[b] / ffeature_smooth_sum
            bfreq_smooth = 1.0 * bfeature_smooth[b] / bfeature_smooth_sum
            score = round(math.log(ffreq_smooth / bfreq_smooth, 2), 3)
            ffreq = round(1.0 * fcount[feature][b] / ffeature_sum, 3)
            bfreq = round(1.0 * bcount[feature][b] / bfeature_sum, 3)
            output += '\t'.join(['', str(b), str(score),\
                                 str(fcount[feature][b]), str(bcount[feature][b]),\
                                 str(ffreq), str(bfreq)]) + '\n'

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
