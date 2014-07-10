# coding: utf-8

import math, random
import mirscanModule


# mirscan
# ------------------------------------------------------------------------------
# args: queryList: a list of miRNA hairpin Candidate objects
#       md: scoring matrix; a dictionary of dictionaries, whose first keys are the names
#           of scoring criteria, second keys are the possible output values for that
#           criterion, and values are numeric scores.  for training, only the first keys
#           are important.
#       train: (optional) boolean; if mirscan is being called for training purposes, this variable
#            will be given a non-False value.
#       starts: (optional) list of dictionaries whose keys are organism names and values are start
#               positions for the miRNA in the corresponding hairpin sequence.  index of a dictionary
#               in the 'starts' list must match the index of the corresponding candidate in the 'queryList'.
# returns: all_data: list of dictionaries, one for each item in queryList.  keys are the names of
#               criteria from 'fdict', values are either the assigned scores (train==False) or
#               the actual values returned by the criteria function calls (train==True).
# ------------------------------------------------------------------------------
# for each candidate, this function prepares a set of data to describe various features of the
# candidate and puts them into a dictionary called 'args'.  that dictionary should contain all
# of the information needed by the various functions that are defined below in 'fdict' (this
# function uses the 'fdict' from its parent environment).  each candidate hairpin or set of
# orthologous hairpins is treated as either an individual miRNA candidate (if 'starts' are
# specified) or as a set of candidates, one for each position along the hairpin at which the
# mature miR could start.  for each candidate position, 'args' is fed into mscore, and a score
# is produced.  the best scoring position for a given hairpin is obtained, along with its score,
# and a dictionary of those criteria scores is returned for each candidate.

def mirscan(queryList,md,train=False,starts=False):

    mirLength = 22  ## user may alter this value, but the variable must remain

#####################################################
## USER SHOULD NOT ALTER THIS BOX
#####################################################
    all_data = []                                   #
    for n in range(len(queryList)):                 #
        args = dict()                               #
        qt = queryList[n]                           #
        args['le'] = mirLength                      #
        # this will guarantee that the org list will end up in the same
        # order for every candidate                 #
        args['orgs'] = qt.organisms()               #
        args['orgs'].sort()                         #
                                                    #
        # get the sequences, folds, and alignments  #
        args['seqs'] = map(lambda org: qt.seq(org), args['orgs'])
        args['ifolds'] = mirscanModule.get_folds(args['seqs'])
        args['al'] = get_alignment(args['seqs'])    #
                                                    #
#####################################################


        ### DEFINE args ENTRIES HERE ###


#####################################################
## USER SHOULD NOT ALTER THIS BOX
#####################################################
                                                    #
        # get the scores for all of the possible (or specified) start positions
        candlist = []                               #
        if starts:                                  #
            args['pos'] = map(lambda org: starts[n][org], args['orgs'])
            args['iside'] = mirscanModule.pick_side(args['ifolds'],args['seqs'],args['pos'],args['le'])
            candlist.append(mirscanModule.mscore(args,md,fdict,train))
        elif train:                                 #
            start_list = mirscanModule.make_start_list(args['seqs'],args['al'],args['le'])
            args['pos'] = random.choice(start_list) #
            args['iside'] = mirscanModule.pick_side(args['ifolds'],args['seqs'],args['pos'],args['le'])
            candlist.append(mirscanModule.mscore(args,md,fdict,train))
        else:                                       #
            start_list = mirscanModule.make_start_list(args['seqs'],args['al'],args['le'])
            for i in start_list:                    #
                args['pos'] = i                     #
                args['iside'] = mirscanModule.pick_side(args['ifolds'],args['seqs'],args['pos'],args['le'])
                candlist.append(mirscanModule.mscore(args,md,fdict))
                                                    #
        # sort the start position options by the totscores and return the highest
        npl = [(x['totscore'],x) for x in candlist] #
        npl.sort()                                  #
        candlist = [val for (key, val) in npl]      #
        bestCand = candlist[-1]                     #
        bestCand['name'] = qt.name                  #
        all_data.append(bestCand)                   #
                                                    #
    return all_data                                 #
#####################################################



# get_alignment
# ------------------------------------------------------------------------------
# args: two sequences
# returns: alignment as list of strings (input sequences + gaps)
# ------------------------------------------------------------------------------
# finds alignment of two sequences
def get_alignment(sa):
    return mirscanModule.global_align(sa[0],sa[1],-8,-3,1,-3,-2,'get traceback')[1:]




# fdict
# this is a dictionary of functions.  each "function" is an instance of class string_ or
# number_feature.  the keys are the names of the criteria.  if the same function is used
# by multiple criteria, then there is a reference for each.
fdict = dict()




### PROTOTYPE fdict ENTRY:
### ------------------------------------------------------------------------------
### First, give the criteria a name (this should be a string, 'nameOfCriteria' below) and
### define it as either a number_feature or a string_feature (these are implemented in
### mirscanModule):
#
# fdict[nameOfCriteria] = mirscanModule.number_feature()
#
### Now, add the two required attributes of the feature: first, a function that measures some aspect
### of the miRNA candidate and returns an appropriate value (must be a string for a string_feature
### and must be a number for a number_feature); and second, a list of all possible values that
### can be returned by the defined function (tip: it may be useful to define the acceptable return
### values first, and then use that list in the function to make sure that the return value is
### acceptable; this will be especially helpful for a number_feature, where generating the return
### value may involve rounding some other quantitative measure.  note also that, if you fail to do
### this, that number_feature has a rounding function that it will apply automatically).
#
# fdict[nameOfCriteria].kl = listOfKeys
# fdict[nameOfCriteria].fx = function that takes 'self' and the argument dictionary that is
#                            defined in the 'mirscan' function as its two arguments.
#
### note that any other set of attributes may be added to these instances, or the
### automatically-assigned attributes may be overwritten.  in the case of adding
### attributes, it is important to respect the namespace of the required attributes;
### in the case of overwriting attributes, it is important to maintain their specified
### requirements.
### ------------------------------------------------------------------------------
