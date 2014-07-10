#!/usr/bin/env python
# coding: utf-8

# This script is used to generate a scoring matrix for evaluating miRNA candidates.
#
# It imports a set of foreground hairpins (a miRNA training set), a set of background
# hairpins (the bulk candidates that the user is trying to filter), and a set of evaluation
# criteria.  It uses the frequencies of each of the possible return values of each criteria
# in the foreground set versus the background set to arrive at a score for that value, returned
# by that feature.
#
# The scoring matrices, along with some additional info about how they were
# generated, are printed out in '.matrix' format.

import math, sys, time, argparse
import mirscanModule as ms

def train(trainfile, criteriafile):
    """
    For each criteria defined in criteriafile, for each possible output value,
    a score is determined which is the base 2 log of the ratio of the frequency
    of that output value in the foreground vs the frequency of that output
    value in the background.  for each criteria, for each output value, the
    foreground and background frequencies are figured such that they include
    pseudocounts that are added according to that feature's 'pseudo' function.

    Output text is in the specified format of a '.matrix' file.
    """
    # the foreground sets
    fqueries,fstarts = ms.parse_train(trainfile, True)

    # the background sets
    bqueries = list()
    for bgf in ms.get_background_files(trainfile):
        bqueries.extend(ms.parse_query(bgf))

    # read the mirscan criteria and build a substitute scoring matrix
    # (mse['bogus_ms']) to serve to the mirscan function when training.
    mse = ms.parse_criteria(criteriafile)
    mse['bogus_ms'] = dict()
    for k in mse['fdict'].keys():
        mse['bogus_ms'][k] = False

    ### for each criteria, get the frequencies of each value
    ### in the foreground and background sets.
    fn = dict()
    bn = dict()
    for k in mse['fdict'].keys():
        fn[k] = dict()
        bn[k] = dict()
        for v in mse['fdict'][k].kl:
            fn[k][str(v)] = 0
            bn[k][str(v)] = 0
    for i in mse['mirscan'](fqueries,mse['bogus_ms'],True,fstarts):
        for j in fn.keys(): fn[j][i[j]]+=1
    for i in mse['mirscan'](bqueries,mse['bogus_ms'],True):
        for j in bn.keys(): bn[j][i[j]]+=1

    ### make header comment lines for the .matrix output with some basic
    ### information about the current training run.
    fcount = sum(fn[mse['fdict'].keys()[0]].values())
    output = '# number = '+str(number)+', fcount = '+str(fcount)+'\n'
    output += '# training file: '+trainfile+'\n'
    output += '# '+time.asctime(time.localtime())+'\n'
    ol = fqueries[0].organisms()
    ol.sort()
    output += '# org keys:\t'+'\t'.join(ol)+'\n'

    ### for each criteria, for each value, generate a score and a line
    ### for the .matrix output that documents that score and the data
    ### that contributed to it.
    for k in mse['fdict'].keys():
        output += k+'\n'
        fcp = mse['fdict'][k].pseudo(fn[k])
        fcpt = sum(fcp.values())
        fct = sum(fn[k].values())
        bcp = mse['fdict'][k].pseudo(bn[k])
        bcpt = sum(bcp.values())
        bct = sum(bn[k].values())
        for v in mse['fdict'][k].kl:
            ff = float(fcp[str(v)])/fcpt
            bf = float(bcp[str(v)])/bcpt
            fraw = round(float(fn[k][str(v)])/fct,3)
            braw = round(float(bn[k][str(v)])/bct,3)
            scr = round(math.log(ff/bf,2),3)
            output += '\t'+str(v)+'\t'+str(scr)+'\t'+str(fn[k][str(v)])+'\t'+str(bn[k][str(v)])+'\t'+str(fraw)+'\t'+str(braw)+'\n'
    return output


parser = argparse.ArgumentParser(description='MiRscan3 Trainer',
            epilog='Paulo Pinto, IEB-WWU, based on:\nhttp://bartellab.wi.mit.edu/softwareDocs/MiRscan3/Introduction.html')

parser.add_argument(dest='trainfile',
                    help='the training file (.train)')
parser.add_argument(dest='criteriafile',
                    help='the mirscan criteria file (.py)')
parser.add_argument('-o', dest='matrixfile', default='stdout',
                    help='the output scoring matrix file (.matrix) (default: stdout)')

args = parser.parse_args()

# check that args corresponding to filenames have the proper extensions.
if args.trainfile.split('.')[-1] != 'train':
    raise ValueError('Training file must be in \'.train\' format.')
if args.criteriafile.split('.')[-1] != 'py':
    raise ValueError('Criteria file must be in \'.py\' format.')
if args.matrixfile.lower() != 'stdout' and args.matrixfile.split('.')[-1] != 'matrix':
    raise ValueError('Output matrix file must be in \'.matrix\' format.')


matrixstring = train(args.trainfile, args.criteriafile)

with (sys.stdout if args.matrixfile.lower() == 'stdout' else open(args.matrixfile, 'w')) as output:
    output.write(matrixstring)
