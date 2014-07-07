import math, random
import mirscanModule

# argv[1] is the query file (.wsl or .fax)
# argv[2] is the scoring matrix file
# argv[3] is the output file

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
        args['orgs'] = qt.orgList()                 #
        args['orgs'].sort()                         #
                                                    #
        # get the sequences, folds, and alignments  #
        args['seqs'] = map(lambda org: qt.seq(org), args['orgs'])
        args['ifolds'] = mirscanModule.get_folds(args['seqs'])
        args['al'] = get_alignment(args['seqs'])    #
                                                    #
#####################################################


        # prepare the basepair dictionaries
        args['bpll'],args['bprl'],args['bpdl'],args['bulgesl'] = [],[],[],[]
        for i in args['ifolds']:
            newbpl,newbpr,newbpd = mirscanModule.make_bp_dicts(i)
            args['bpll'].append(newbpl)
            args['bprl'].append(newbpr)
            args['bpdl'].append(newbpd)
            args['bulgesl'].append(method5_pre_help(i,newbpd))

        # get the complexity measurement for the hairpin
        args['compl'] = map(lambda a: max([complexity_pre_help(a,len(a)/5), \
                                   complexity_pre_help(mirscanModule.reverse_string(a), \
                                                       len(a)/5)]),args['seqs'])


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
        bestCand['name'] = qt.name()                #
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




# method5_alt_D 'bulge_sym_D'
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# help args: seq: a sequence (or one sequence from an alignment)
#            blg: a list of bulges as generated by method5_pre_help for the 'al' entry of interest
#            pos: integer starting position of the miR in the 'al' entry of interest
#            le: integer imposed miR length
# help returns: integer number of unsymmetrically-bulged bases in the miR
# returns: arithmatic mean of help returned values over all candidate ortholog seqs
# ------------------------------------------------------------------------------
# for each bulge, counts number of unsymmetrically-bulged bases and sums them over all bulges
# uses individual structures
fdict['bulge_sym_D'] = mirscanModule.number_feature()
fdict['bulge_sym_D'].kl = [0,1,2,3,4,5,6,7,8,9]
def method5_alt_D(self,ar):
    def fx5_help(seq,blg,pos,le):
        a = b = asym = 0
        while a<le-2:
            if seq[pos+b]!='-': a+=1
            b+=1
        for i in blg:
            if i.b.bef+1>=pos and i.b.aft<=pos+b:
                ca = i.b.aft - i.b.bef - seq[i.b.bef+1:i.b.aft].count('-')
                cb = i.o.aft - i.o.bef - seq[i.o.bef+1:i.o.aft].count('-')
                asym += abs(ca-cb)
        return asym
    return float(sum(map(lambda r: fx5_help(ar['seqs'][r],ar['bulgesl'][r],ar['pos'][r],ar['le']), range(len(ar['al'])))))/len(ar['al'])
fdict['bulge_sym_D'].fx = method5_alt_D


# method5_pre_help
# ------------------------------------------------------------------------------
# args: fold: string of the bracket-notation nucleic acid 2ndary structure
#       bpd: base pair dictionary, as generated by make_bp_dicts
# returns: bulges: a list of bulges as described below
# ------------------------------------------------------------------------------
# creates a list of bulges for the hairpin to be passed to function fx5 whenever it is
# called for that hairpin.  each element in the list is a bulge; each bulge is an instance
# of class bpoint (defined here), which defines the bulge as internal state "b" and the
# opposite side of the bulge as internal state "o".  internal states "b" and "o" are each
# instances of class b2 (defined here), which has internal state "bef" (position of the
# nucleotide before the opening of the bulge) and internal state "aft" (position of the
# nucleotide which is paired to close the bulge).  unpaired positions at the beginning and
# end of the sequence are not counted as bulges.  all bulges appear twice in the list (such
# that the internal states "b" and "o" are switched).  loops appear once, with their "bef"
# and "aft" values stored in the internal state "b".  internal state "o" carries the value
# of 0 as its "bef" and "aft" internal states in the case of loops.
def method5_pre_help(fold,bpd):
    foldlength = len(fold)
    bulgesl = []
    bulgel = []
    for n in range(foldlength):
        bulgel.append(n)
        if fold[n]!='.':
            bulgesl.append(bulgel)
            bulgel = []
    bulgel.append(foldlength)
    bulgesl.append(bulgel)
    class bpoint:
        b = False # bulge
        o = False # other
    class b2:
        bef = False
        aft = False
    bulges = []
    bulges_first = []
    bulges_help = []
    for i in fold: bulges_help.append(bpoint())
    for i in bulgesl:
        new = b2()
        if min(i)>0: new.bef = min(i)-1
        if max(i)<foldlength: new.aft = max(i)
        if new.bef and new.aft:
            bulges_first.append(bpoint())
            bulges_first[len(bulges_first)-1].b = new
            for j in range(min(i),max(i)): bulges_help[j].b = new
    loop_opp = b2()
    loop_opp.bef = loop_opp.aft = 0
    for i in bulges_first:
        if i.b.bef!=bpd[i.b.aft]:
            if bpd[i.b.bef]>i.b.bef and bulges_help[bpd[i.b.bef]-1].b: i.o = bulges_help[bpd[i.b.bef]-1].b
            elif bpd[i.b.bef]<i.b.bef and bulges_help[bpd[i.b.bef]-1].b: i.o = bulges_help[bpd[i.b.bef]-1].b
            else:
                i.o = b2()
                i.o.bef = bpd[i.b.aft]
                i.o.aft = bpd[i.b.bef]
        else: i.o = loop_opp
        if i.o and i.b.aft-i.b.bef+i.o.aft-i.o.bef>2: bulges.append(i)
    return bulges


# method6_G
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
#       give_edges: controls the returned value (optional; see below)
# help args: fold: string of the bracket-notation nucleic acid 2ndary structure
#            seq: a sequence (or one sequence from an alignment)
#            bpd: base pair dictionary, as generated by make_bp_dicts
#            pos: integer starting position of the miR in the 'al' entry of interest
#            le: integer imposed miR length
#            side: 'left' if the miR is on the left (5') side of the hairpin;
#                  'right' if the miR is on the right (3') side of the hairpin;
#                   as decided by pick_side
#            give_edges: controls the returned value (optional; see below)
# help returns: if give_edges is False (not given), returns the length of the loop in nucleotides
#               (loop defined as the stretch between the miR and miR*).  otherwise, returns a two-item
#               list with the start and end coordinates of the loop
# returns: arithmatic mean of help returned values (loop lengths) over all candidate ortholog seqs
#          if give_edges is not False, then it returns a list of two-item lists, where the
#          outer list items are folds from args['ifolds'] and the inner list has the start and
#          end coordinates for that fold's loop (loop defined as the stretch between the miR and miR*)
# ------------------------------------------------------------------------------
# measures the loop size (loop defined as the stretch between the miR and miR*)
# using ar['iside'] and individual structures
# an implementation for counting loop length which is more ornate
# using ar['iside'] and individual structures
# altered from old 6_E to no longer need to look for '-'
def method6_G(self,ar,give_edges=False):
    def fx6_help(fold,seq,bpd,pos,le,side,give_edges=False):
        brackets = { '(':1, '.':0, ')':-1, 'X':0 }
        p = pos+le
        beginning = end = False
        if side=='left':
            x = 0
            while pos+x<p and fold[pos+x]!='(': x+=1
            if pos+x<p:
                beginning = p
                end = bpd[pos+x]-(le-1-x-2)
        elif side=='right':
            x = 1
            while p-x>=pos and fold[p-x]!=')': x+=1
            if p-x>=pos:
                beginning = bpd[p-x]+(le-x+3)
                end = pos
        if not(beginning):
            def edge_helper(pos,p,el,er,seq,le,side):
                a = pos - el - 2
                b = er - p
                if side=='right':
                    leng = a - b - le
                    return pos-(leng-2),pos
                else:
                    leng = b - a - le
                    return p,p+(leng+2)
            def left_helper(pos,p,fold,seq,le,side):
                x = 0
                while fold[p+x]!='(': x+=1
                if side=='left': return p,bpd[p+x]-(x+3)
                else: return edge_helper(pos,p,0,p+x,seq,le,side)
            def right_helper(pos,p,fold,seq,le,side):
                x = 1
                while fold[pos-x]!=')': x+=1
                if side=='right':
                    if x-3<0: return bpd[pos-x]+(1-x+3),pos
                    else: return bpd[pos-x]-(x-3),pos
                else: return edge_helper(pos,p,pos-x+1,len(seq),seq,le,side)
            sides = [filter(lambda a: a!='.', 'X'+fold[:pos])[-1],filter(lambda a: a!='.', fold[p:]+'X')[0]]
            if sides==['(',')']:
                el = pos-1
                while fold[el]=='.': el-=1
                er = p
                while fold[er]=='.': er+=1
                beginning,end = edge_helper(pos,p,el,er,seq,le,side)
            elif sides==['(','(']:
                el = pos-1
                while fold[el]=='.': el-=1
                beginning,end = edge_helper(pos,p,el,bpd[el],seq,le,side)
            elif sides==[')',')']:
                er = p
                while fold[er]=='.': er+=1
                beginning,end = edge_helper(pos,p,bpd[er],er,seq,le,side)
            elif sides==['X','(']: beginning,end = left_helper(pos,p,fold,seq,le,side)
            elif sides==[')','X']: beginning,end = right_helper(pos,p,fold,seq,le,side)
            elif sides==[')','(']:
                if fold[:pos].count('(')+fold[:pos].count(')') < fold[p:].count('(')+fold[p:].count(')'):
                    beginning,end = left_helper(pos,p,fold,seq,le,side)
                else: beginning,end = right_helper(pos,p,fold,seq,le,side)
            elif sides==['X','X']: beginning,end = edge_helper(0,len(fold),0,len(fold),seq,le,side)
        if give_edges: return [beginning,end]
        else: return end - beginning
    if give_edges: return map(lambda r: fx6_help(ar['ifolds'][r],ar['seqs'][r],ar['bpdl'][r],ar['pos'][r],ar['le'],ar['iside'],1), range(len(ar['al'])))
    else: return float(sum(map(lambda r: fx6_help(ar['ifolds'][r],ar['seqs'][r],ar['bpdl'][r],ar['pos'][r],ar['le'],ar['iside']), range(len(ar['al'])))))/len(ar['al'])
fdict['loop_dis_G'] = mirscanModule.number_feature()
fdict['loop_dis_G'].kl = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]
fdict['loop_dis_G'].fx = method6_G


# method7
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
#            uses ar['seqs'] and ar['pos']
# returns: a one-letter string
# ------------------------------------------------------------------------------
# returns the one-letter string of the nucleotide at the position number of the sequence
# string provided.  'N' appears in genome sequence files and is also included in the nucs
# list.  if another character appears, it is printed and returned (an exception will result).
def method7(self,ar):
    if self.kl.count(ar['seqs'][self.seq][ar['pos'][self.seq]+self.pos]) == 0: print ar['seqs'][ar['pos'][self.seq]]
    return ar['seqs'][self.seq][ar['pos'][self.seq]+self.pos]

for n in [0,8]:
    for nn in range(2):
        name = 'nuc'+str(n+1)+'_s'+str(nn+1)
        fdict[name] = mirscanModule.string_feature()
        fdict[name].pos = n
        fdict[name].seq = nn
        fdict[name].kl = ['A','T','C','G','N']
        fdict[name].fx = method7


# complexity
# from lempel & ziv (IEE proc inf theory IT-22: 75-81; 1976) and
# denominator is from gusev et al (biosystems 30:200-231; 1993)
# figured initially; method just averages the initial results of all the input seqs.
#
# method_complexity
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
#            uses ar['compl']
# returns: the arthmatic mean of the values in the ar['compl_A12'] list
# ------------------------------------------------------------------------------
def method_complexity(self,ar): return sum(ar['compl'])/len(ar['compl'])
fdict['compl'] = mirscanModule.number_feature()
fdict['compl'].kl = [-2,-1,0,1,2,3,4,5]
fdict['compl'].fx = method_complexity
#
# complexity_pre_help
# ------------------------------------------------------------------------------
# args: seq: a string (a sequence of letters)
#       level_max: an integer representing the longest repeated word that will be
#            sought
# returns: a normalization of the longest repeated word's length to the length of the
#               overall sequence.  for sequence which is actually random, this value
#               remains approximately constant with length of the sequence
# ------------------------------------------------------------------------------
def complexity_pre_help(seq,level_max):
    tree = dict()
    sl = len(seq)
    pos = c = 0
    steps = []
    words = []
    while pos<sl:
        step = find_in_tree(tree,1,pos,seq,level_max,sl)
        words.append(seq[pos:pos+step])
        pos += step
        steps.append(step)
        for n in range(1,level_max+2):
            if pos-n>=0:
                find_in_tree(tree,1,pos-n,seq,n,sl,1)
        c+=1
    longest = max(steps)
    return math.log(float(4**longest)/(len(seq)-longest*2)**2,4)
#
# find_in_tree
# ------------------------------------------------------------------------------
# args: tree: a series of nested dictionaries representing all of the sequences
#            that have been seen prior in the sequence, or a branch from that tree
#            keys are nucleotide letters; values are dictionaries with key:value pairs
#            for letters that appear after that letter (or series of letters) somewhere
#            in the sequence
#       level: integer; how many steps into the nested dictionary have been taken
#       pos: integer; position in seq that is currently being compared to the tree
#       seq: string; sequence being examined
#       level_max: integer; maximum number of steps into the tree that will be taken
#            (this is the same as the max repeated word length)
#       sl: integer; length of 'seq'
#       build: optional; if not false, updates (mutates) 'tree' to include the whole
#            sequence up until the new position
# returns: integer corresponding to the number of steps that have been taken into the
#               tree (ie the number of steps forward that 'pos' should be advanced in
#               'seq' to get to the new starting position
# ------------------------------------------------------------------------------
def find_in_tree(tree,level,pos,seq,level_max,sl,build=False):
    if pos==sl: return level - 1
    elif tree.has_key(seq[pos]):
        if level==level_max: return level
        else: return find_in_tree(tree[seq[pos]],level+1,pos+1,seq,level_max,sl,build)
    else:
        if build:
            tree[seq[pos]] = dict()
            if level!=level_max: return find_in_tree(tree[seq[pos]],level+1,pos+1,seq,level_max,sl,build)
        if level==level_max: return level-1
        else: return max([level-1,1])


# bp_matrix_A
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# help args: fold: string of the bracket-notation nucleic acid 2ndary structure
#            pos: integer starting position of the miR in the 'al' entry of interest
#            side: 'left' if the miR is on the left (5') side of the hairpin;
#                  'right' if the miR is on the right (3') side of the hairpin;
#                   as decided by pick_side
# help returns: 1 or -1 (depending on direction of the pair) or 0 (for no pair)
# returns: string 'paired' if position is paired in all structures; 'unpaired' otherwise
# ------------------------------------------------------------------------------
# these functions use individual structure output and evaluate the probabilities of base
# pairs at specific sequence positions.  if the position is paired, the helper function
# returns an integer '1' (left bracket) or '-1' (right bracket); if not, an integer '0'.
# the function retrurns the absolute value of the average of those values, i.e. the fraction
# of sequences with that position paired.  for each instance, the position # relative to the
# miR start is the 'pos' internal state.
def method8_A(self,ar):
    def fx1_help(fold,pos,side):
        fd = { '(':1, ')':-1, '.':0 }
        if pos+self.pos<0 or pos+self.pos>=len(fold): return 0
        elif side=='left': return fd[fold[pos+self.pos]]
        elif side=='right': return -fd[fold[pos+self.pos]]
    if float(sum(map(lambda r: fx1_help(ar['ifolds'][r],ar['pos'][r],ar['iside']), range(len(ar['al'])))))/len(ar['al'])==1: return 'paired'
    else: return 'unpaired'

for n in range(0,20):
    if n >= 0:
        entry_name = 'bp_matrix_A_n'+str(n+1)
        fdict[entry_name] = mirscanModule.string_feature()
        fdict[entry_name].kl = ['paired','unpaired']
        fdict[entry_name].pos = n
        fdict[entry_name].fx = method8_A
    else:
        entry_name = 'bp_matrix_A_n'+str(n)
        fdict[entry_name] = mirscanModule.string_feature()
        fdict[entry_name].kl = ['paired','unpaired']
        fdict[entry_name].pos = n
        fdict[entry_name].fx = method8_A


# bp_matrix_C
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# help args: fold: string of the bracket-notation nucleic acid 2ndary structure
#            pos: integer starting position of the miR in the 'al' entry of interest
#            side: 'left' if the miR is on the left (5') side of the hairpin;
#                  'right' if the miR is on the right (3') side of the hairpin;
#                   as decided by pick_side
# help returns: 1 or -1 (depending on direction of the pair) or 0 (for no pair)
# returns: string 'paired' if position is paired in all structures; 'unpaired' otherwise
# ------------------------------------------------------------------------------
# these functions use individual structure output and evaluate the probabilities of base
# pairs at specific sequence positions.  if the position is paired, the helper function
# returns an integer '1' (left bracket) or '-1' (right bracket); if not, an integer '0'.
# the function retrurns the absolute value of the average of those values, i.e. the fraction
# of sequences with that position paired.  for each instance, the position # relative to the
# miR start is the 'pos' internal state.
# uses individual structures
# changes orientation of sequences on the right side of the hairpin so that we can
# see pairing probability as a function of position within the miR-miR* duplex
def method8_C(self,ar):
    def fx1_help(fold,pos,le,side):
        fd = { '(':1, ')':-1, '.':0 }
        if side=='left':
            if pos+self.pos<0 or pos+self.pos>=len(fold): return 0
            else: return fd[fold[pos+self.pos]]
        elif side=='right':
            if pos+le-2-self.pos>=len(fold) or pos+le-2-self.pos<0: return 0
            else: return -fd[fold[pos+le-2-self.pos]]
    if float(sum(map(lambda r: fx1_help(ar['ifolds'][r],ar['pos'][r],ar['le'],ar['iside']), range(len(ar['al'])))))/len(ar['al'])==1: return 'paired'
    else: return 'unpaired'

for n in [-8,-7,-6,-5,-4,-3,-2,-1,20,21]:
    if n >= 0:
        entry_name = 'bp_matrix_C_n'+str(n+1)
        fdict[entry_name] = mirscanModule.string_feature()
        fdict[entry_name].kl = ['paired','unpaired']
        fdict[entry_name].pos = n
        fdict[entry_name].fx = method8_C
    else:
        entry_name = 'bp_matrix_C_n'+str(n)
        fdict[entry_name] = mirscanModule.string_feature()
        fdict[entry_name].kl = ['paired','unpaired']
        fdict[entry_name].pos = n
        fdict[entry_name].fx = method8_C


# method_con: conservtion matrix
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# returns: all_match: string 'con' if the nuc identity at the positin given by
#               self.pos is the same in all sequences; otherwise, returns 'non'
# ------------------------------------------------------------------------------
def method_con(self,ar,more_args=False):
    seqs = ar['seqs']
    if more_args:
        pos = map(lambda n: n+more_args['start'], ar['pos'])
        all_match = 0
        for n in range(more_args['length']):
            add_this = 1
            for n in range(1,len(seqs)):
                if seqs[0][pos[0]]!=seqs[n][pos[n]]: add_this = 0
            all_match += add_this
        return all_match
    else:
        pos = map(lambda n: n+self.pos, ar['pos'])
        all_match = 'con'
        for n in range(1,len(seqs)):
            if seqs[0][pos[0]]!=seqs[n][pos[n]]: all_match = 'non'
        return all_match

for n in [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21]:
    fdict['con'+str(n+1)] = mirscanModule.string_feature()
    fdict['con'+str(n+1)].kl = ['con','non']
    fdict['con'+str(n+1)].pos = n
    fdict['con'+str(n+1)].fx = method_con


# loop_vs_mir
# uses method9_B (below) to figure out the fraction of the loop that is conserved;
# uses method_con with more_args to get the fraction of the miR that is conserved;
# returns the difference between those two fractions
#
# method343
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# returns: float; difference between fraction of miR conserved and fraction of loop
#               conserved; positive if miR is more conserved
# ------------------------------------------------------------------------------
# compares conservation of miR with conservation of loop; uses method_con
def method343(self,ar):
    more_args = {'length':22, 'start':0}
    mircon = float(method_con(self,ar,more_args))/more_args['length']
    loopcon = method9_B(self,ar)
    return mircon - loopcon
#
# method9_B
# ------------------------------------------------------------------------------
# args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# returns: float; # conserved nucs / # total nucs in loop
# ------------------------------------------------------------------------------
# looks at loop conservation; uses method6_G
def method9_B(self,ar):
    loop_list = method6_G(self,ar,1)
    seq_list = map(lambda r: ar['seqs'][r][loop_list[r][0]:loop_list[r][1]+1],range(len(loop_list)))
    if max(map(len,seq_list))==0: return 0
    else:
        al = get_alignment(seq_list)
        match=0
        for n in range(len(al[0])):
            if al[0][n]==al[1][n]: match+=1
        return float(match)/max(map(len,al))

fdict['loop_vs_mir'] = mirscanModule.number_feature()
fdict['loop_vs_mir'].kl = map(lambda a: .05*a, range(-20,21))
fdict['loop_vs_mir'].fx = method343
