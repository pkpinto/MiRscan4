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
# The scoring matrices, along with some additional info about
# how they were generated, are printed out in '.matrix' format.  'train'
# is the master function, and is called at the bottom of the script.
#
# argv[1] is the .train training file with the foreground set and the background files
# argv[2] is the mirscan script that will be used for training
# argv[3] is the number of randomly-selected background seqs which will contribute to the score
#         (if set to 0, all background sequences will contribute)
# argv[4] is a .matrix file with the scoring matrix (optional; printed to stdout if left blank) 


import math, sys, random, time, mirscanModule


trainFile = sys.argv[1]
criteriaFile = sys.argv[2]
number = int(sys.argv[3])
matrixFile = sys.argv[4]

# check that args corresponding to filenames have the proper extensions.
if trainFile.split('.')[-1]!='train': raise ValueError('training file must be ".train" format')
if criteriaFile.split('.')[-1]!='py': raise ValueError('criteria file must be formatted for python (".py")')
if matrixFile.lower()!='stdout' and matrixFile.split('.')[-1]!='matrix':
    raise ValueError('outfile must be ".matrix".')





# train
# ------------------------------------------------------------------------------
# args: trainfile: name of a .train format training file
#       msfile: name of a .py format mirscan criteria file
#       number: the number of background hairpins to use for training
#               (0 -> use all available background hairpins)
# returns: string with the text of the score matrix ('.matrix'-formatted)
# ------------------------------------------------------------------------------
# for each criteria defined in msfile, for each possible output value, a score is determined
# which is the base 2 log of the ratio of the frequency of that output value in the foreground
# vs the frequency of that output value in the background.  for each criteria, for each output
# value, the foreground and background frequencies are figured such that they include pseudocounts
# that are added according to that feature's 'pseudo' function.  output text is in the specified
# format of a '.matrix' file.

def train(trainfile,msfile,number):
    ### gather the mirscan criteria and build a substitute scoring matrix
    ### (mse['bogus_ms']) to serve to the 'mirscan' function when training.
    mse = mirscan_eval(msfile)
    mse['bogus_ms'] = dict()
    for k in mse['fdict'].keys():
        mse['bogus_ms'][k] = False

    ### get the foreground and background sets
    fqueries,fstarts = mirscanModule.parseTrain(trainfile,True)
    bqueries = getBackground(trainfile,number)

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
    ol = fqueries[0].orgList()
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



def getBackground(trainfile,number):
    backgroundFiles = mirscanModule.getBackgroundFiles(trainfile)
    bqueries = []
    if number==0:
        for bgf in backgroundFiles: bqueries.extend(mirscanModule.get_queries(bgf))
    else:
        bgfToSize = dict()
        for bgf in backgroundFiles: bgfToSize[bgf] = len(mirscanModule.get_queries(bgf))
        if sum(bgfToSize.values()) < number:
            raise ValueError("not enough queries ("+str(sum(bgfToSize.values()))+") for background of "+str(number))
        else:
            entries = range(sum(bgfToSize.values()))
            random.shuffle(entries)
            useEntries = entries[:number]
            useEntries.sort()
            bottomNum = 0
            for bgf in backgroundFiles:
                theseEntries = filter(lambda n: bottomNum <= n < bottomNum + bgfToSize[bgf], useEntries)
                theseEntries = map(lambda n: n - bottomNum, theseEntries)
                if len(theseEntries) > 0:
                    theseQueries = mirscanModule.get_queries(bgf)
                    bqueries.extend(map(lambda n: theseQueries[n], theseEntries))
                bottomNum += bgfToSize[bgf]
    if number!=0 and number!=len(bqueries): raise ValueError("number "+str(number)+" didn't match num queries: "+str(len(bqueries)))
    return bqueries




# mirscan_eval
# ------------------------------------------------------------------------------
# args: name of a mirscan criteria file (.py format)
# returns: an environment (dictionary) generated by evaluating the input file
# requires: any modules that the input file requires must be accessible from
# the working directory from which THIS training script is being run.
# ------------------------------------------------------------------------------
def mirscan_eval(filename):
    openfile = open(filename)
    code = openfile.read()
    openfile.close()
    envi = dict()
    exec code in envi
    return envi





scorestring = train(trainFile,criteriaFile,number)

if matrixFile.lower()=='stdout':
    sys.stdout.write(scorestring)
else:
    aa = open(matrixFile,'w')
    aa.write(scorestring)
    aa.close()
