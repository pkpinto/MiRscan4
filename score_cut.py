#!/usr/bin/env python
# coding: utf-8

# this script collects the scores earned by the foreground mirs and finds
# the mean and standard dev. of the foreground score distribution.  it decides
# where to cut the background.  this point will be the score which is one
# standard deviation below the lowest score from the foreground set.  it then
# filters through all of the .scr files in the directory and generates new
# .wsl files in a new directory with those sequences that pass the threshold
# score.

# argv[1] is the foreground .scr file
# argv[2] is the background .scr file
# argv[3] (optional) is the .fax background file that was being scored.
# argv[4] (optional) is the new .fax file to which the passing entries will
#         be written.  if 'stdout', .fax-format text will be written to
#         standard out.

import sys, math, mirscanModule


foregroundScores = sys.argv[1]
backgroundScores = sys.argv[2]
if foregroundScores.split('.')[-1]!='scr':
    raise ValueError("foreground scores must be '.scr' format, not "+foregroundScores)
if backgroundScores.split('.')[-1]!='scr':
    raise ValueError("background scores must be '.scr' format, not "+backgroundScores)

if len(sys.argv)>4:
    makeNewFile = True
    sourceFile = sys.argv[3]
    outfileName = sys.argv[4]
    if sourceFile.split('.')[-1]!='fax' and sourceFile.split('.')!='train':
        raise ValueError("background source file must be '.fax' format, not "+sourceFile)
    if outfileName.split('.')[-1]!='fax':
        raise ValueError("filtered output file must be '.fax' format, not "+outfileName)
else:
    makeNewFile = False




# this list can be modified to change the behavior of score_cut.py, but here,
# it is set to use the pre-computed (by mirscan) sum of all of the individual
# feature scores, referred to as 'totscore' in mirscan's output.
keys = ['totscore']


# filtered_score
# -----------------------------------------------------------------
# args: mir: a dictionary of score features (aka scores for features) for
#            a particular miR
#       kl: list of strings which correspond to 'mir' dictionary keys;
#            these are the keys whose values will be summed to generate
#            the final score
# returns: new: sum of the scores for features listed in 'kl'
# -----------------------------------------------------------------
# used to get the scores for individual mirs
def filtered_score(mir,kl):
    new = 0
    for i in kl: new+=mir[i]
    return new



# find_stats
# -----------------------------------------------------------------
# args: mirs: dictionary of dictionaries; output 'c' from get_scores
#       kl: list of strings which correspond to 'mir' dictionary keys;
#            these are the keys whose values will be summed to generate
#            the final score
# returns: mean: arithmatic mean of the scores for all of the mirs in
#               'mirs' arg as determined by filtered_score given 'kl' arg
#          stdev: standard deviation of scores from 'mean'
#          min(a): smallest score for all 'mirs' entries
# -----------------------------------------------------------------
# gets mean score, stdev of scores, and min score for all of the miR
# candidates' filtered scores
def find_stats(mirs,kl):
    a = map(lambda m: filtered_score(m,kl), mirs)
    mean = sum(a)/len(a)
    dev = 0
    for n in a:
        diff = n-mean
        dev += diff*diff
    stdev = math.sqrt(dev/len(a))
    return mean,stdev,min(a)



# get_cands_above_x
# -----------------------------------------------------------------
# args: x: a number
#       mirs: dictionary of dictionaries; output 'c' from get_scores
#       kl: list of strings which correspond to 'mir' dictionary keys;
#            these are the keys whose values will be summed to generate
#            the final score
# returns: new: list of keys from arg 'mirs'
# -----------------------------------------------------------------
# gets a list of all the candidates whose filtered scores are above arg 'x'
def get_cands_above_x(x,mirs,kl):
    new = []
    for m in mirs:
        if filtered_score(m,kl) > x: new.append(m)
    return new





# below, the score distribution for the foreground set is analyzed to
# select a score threshold ('cut'), and the candidates are filtered for
# having scores above the threshold ('babove' are those candidtates' score
# dictionaries).
fs = mirscanModule.get_scores(foregroundScores,keys)
bs = mirscanModule.get_scores(backgroundScores,keys)
fmean,fstdev,fmin = find_stats(fs,keys)
cut = fmin - fstdev/2
babove = get_cands_above_x(cut,bs,keys)


print '\t'.join(['fore','mean='+str(fmean),'stdev='+str(fstdev)])
print 'cut = ',cut
print 'candidates above minimum: ',len(babove),' of ',len(bs)


# if the appropriate arguments are provided, the passing candidates will
# be put into a new .fax file.
if makeNewFile:
    candList = mirscanModule.get_queries(sourceFile)
    passedNames = dict()
    for b in babove: passedNames[b['name']] = None
    passedList = filter(lambda c: passedNames.has_key(c.name()), candList)
    if len(passedList)!=len(babove):
        raise ValueError("lengths of scored ("+str(len(babove))+") and queried ("+\
                         str(len(passedList))+") candidates don't match after filtering.")
    mirscanModule.writeFax(passedList,outfileName)
