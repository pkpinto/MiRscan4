#!/usr/bin/env python
# coding: utf-8

import math, sys, time, argparse
import mirscanIO as msio

def train(forefile, backfile, criteriafile):
    """
    For each criteria defined in criteriafile, for each possible output value,
    a score is determined which is the base 2 log of the ratio of the frequency
    of that output value in the foreground vs the frequency of that output
    value in the background. For each criteria, for each output value, the
    foreground and background frequencies are figured such that they include
    pseudocounts that are added according to that feature's 'pseudo' function.

    Output text is in the specified format of a '.matrix' file.
    """
    # the foreground sets
    fqueries = msio.parse_fam(forefile)
    fstarts = list()
    for c in fqueries:
        __starts = dict()
        for org in c.organisms():
            __starts[org] = c.hairpin(org).find(c.mature(org))
        fstarts.append(__starts)

    # the background sets
    bqueries = msio.parse_query(backfile)

    # read the mirscan criteria and build a substitute scoring matrix
    # (mse['bogus_ms']) to serve to the mirscan function when training
    mse = msio.parse_criteria(criteriafile)
    mse['bogus_ms'] = dict()
    for k in mse['fdict'].keys():
        mse['bogus_ms'][k] = False

    # for each criteria, get the frequencies of each value in the foreground
    # and background sets
    fn = dict()
    bn = dict()
    for k in mse['fdict'].keys():
        fn[k] = dict()
        bn[k] = dict()
        for v in mse['fdict'][k].kl:
            fn[k][str(v)] = 0
            bn[k][str(v)] = 0
    for i in mse['mirscan'](fqueries, mse['bogus_ms'], True, fstarts):
        for j in fn.keys(): fn[j][i[j]]+=1
    for i in mse['mirscan'](bqueries, mse['bogus_ms'], True):
        for j in bn.keys(): bn[j][i[j]]+=1

    # make header comment lines for the .matrix output with some basic
    # information about the current training run
    fcount = sum(fn[mse['fdict'].keys()[0]].values())
    output = '# number = 0, fcount = '+str(fcount)+'\n'
    output += '# training file: '+forefile+'\n'
    output += '# '+time.asctime(time.localtime())+'\n'
    ol = fqueries[0].organisms()
    ol.sort()
    output += '# org keys:\t'+'\t'.join(ol)+'\n'

    # for each criteria, for each value, generate a score and a line for the
    # .matrix output that documents that score and the data that originated it
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


matrixstring = train(args.forefile, args.backfile, args.criteriafile)

with (sys.stdout if args.matrixfile.lower() == 'stdout' else open(args.matrixfile, 'w')) as output:
    output.write(matrixstring)
