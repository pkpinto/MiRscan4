# coding: utf-8

import string, os

# Hairpin criteria objects (features)
# ------------------------------------------------------------------------------
# these objects are used more like dictionaries, simply to cluster together attributes
# that are related to the evaluation of a particular hairpin feature.  the user may
# add additional attributes to a feature object as needed; the requirement is that
# each object, in addition to having the attributes implemented below, also be
# given an attribute called 'fx' and another called 'kl' (described below).
#
# all f(x) methods for obtaining values associated with a scoring criteria will
# be the "fx" method for an object of the criteria_measure class, which will inherit
# from an instance of either the StringFeature or NumericalFeature class.  each of those
# classes has an internal state "type", which returns either 'string' or 'number' as
# a string, "kl" which is a list of possible key values from the scoring tables, and
# "kv", which is a method for taking a value returned by a "fx" method and turning it into
# an appropriate key value.  note that numerical values in kl must be in ascending order.
# there is also a function "pseudo" for each of the types of features which takes a string
# of counts and adds pseudocounts.
#
# StringFeature
# ------------------------------------------------------------------------------
# attributes: fx: function for evaluating the feature in question (initially, boolean False);
#                  must take self and args dictionary as only arguments and return a string
#             type: string labelling the instance as a 'string'-type feature
#             kl: list of all possible return values that are accepted for that feature
#             kv: function for binning string returned by fx.  for strings, the arg value x
#                  is just returned
#             ex: function
#             pseudo: function
# ------------------------------------------------------------------------------
class StringFeature:
    type = 'string'
    def kv(self,x):
        return x
    def ex(self,*args):
        return self.kv(self.fx(*args))
    def pseudo(self,cd):
        new = dict()
        for k in cd.keys():
            new[k] = cd[k]+1
        return new

# NumericalFeature
# ------------------------------------------------------------------------------
# attributes: fx: function for evaluating the feature in question (initially, boolean False);
#                  must take self and args dictionary as only arguments and return a number
#             type: string labelling the instance as a 'number'-type feature
#             kl: list of all possible return values (in ascending numerical order)
#                  that are accepted for that feature
#             kv: function for binning number returned by fx.  for numbers, the arg value x
#                  is placed in the kl bin where the kl value is less than or equal to x and
#                  where the next bin's value (if there is one) is greater than x
#             ex: function
#             pseudo: function
# ------------------------------------------------------------------------------
class NumericalFeature:
    type = 'number'
    def kv(self,x):
        p = 1
        while p<len(self.kl) and self.kl[p]<=x: p+=1
        return str(self.kl[p-1])
    def ex(self,*args):
        return self.kv(self.fx(*args))
    def pseudo(self,cd,iter=0):
        kernel = [.75,.125]
        new = dict()
        for k in cd.keys(): new[k]=0
        if iter==2:
            for k in cd.keys():
                new[k] = cd[k]+1
            return new
        else:
            for n in range(len(self.kl)):
                new[str(self.kl[n])]+=kernel[0]*cd[str(self.kl[n])]
                for m in range(1,len(kernel)):
                    for s in [-1,1]:
                        if n+s*m>=0 and n+s*m<len(self.kl): new[str(self.kl[n+s*m])]+=kernel[m]*cd[str(self.kl[n])]
            return self.pseudo(new,iter+1)



# get_folds
# ------------------------------------------------------------------------------
# args: sl: a list of sequences
# returns: ind_folds: a list of bracket-notation folds corresponding to the mfe
#               fold for the sequence of the same index in the sl list
# requires: RNAfold (vienna package) is installed and executable at the command
#               line by typing 'RNAfold' to begin an interactive session.
# ------------------------------------------------------------------------------
# uses RNAfold (zuker algorithm) to generate an mfe fold for each of the sequences
# in sl
def get_folds(sl):
    ind_folds = []
    if len(sl)>0:
        ut_table = string.maketrans('U','T')
        all_s = '\n'.join(sl)
        fi,fo = os.popen2('RNAfold --noPS')
        fi.write(all_s)
        fi.close()
        b = fo.read().split('\n')
        fo.close()
        b_index = 0
        for s in sl:
            if b[b_index].translate(ut_table)==s:
                c = b[b_index+1].split(' ')
                b_index+=2
            else:
                c = ['.'*len(s)]
            ind_folds.append(c[0])
    return ind_folds



# make_start_list
# ------------------------------------------------------------------------------
# args: seq_list: a list of the query sequences with no gaps ('-' characters)
#       al: a list of the query sequences as they would appear in an alignment
#            with one another (ie with gaps, all strings should be the same length)
#       length: integer length of a miR
# returns: pos_list: a list of lists.  each element in the outer list corresponds
#               to a set of start positions for a miR.  the first n elements of the
#               the inner list gives the starting position in each of the n sequences
#               of seq_list, in the same order as the sequences they refer to.  the
#               n+1 element gives the starting position in the alignment.
# ------------------------------------------------------------------------------
# generates a list of lists, where the elements of the list are possible start positions for the mir
# in each of the sequences; start positions in the returned list are in the order of the sequences
# within the input sequence list.
def make_start_list(seq_list,al,length):
    sacl = [True]
    pos_list = []
    for n in range(len(seq_list)):
        sacp = dict()
        rp = ap = 0
        while rp<len(seq_list[n])-length+1:
            sacp[ap]=rp
            if al[n][ap]!='-': rp+=1
            ap+=1
        sacl.append(sacp)
    for n in range(len(al[0])):
        if reduce(lambda a,b: a and b.has_key(n), sacl):
            new_entry = map(lambda a: a[n], sacl[1:])
            new_entry.append(n)
            pos_list.append(new_entry)
    while pos_list[1].count(0)>0:
        pos_list=pos_list[1:]
    return pos_list



# pick_side
# ------------------------------------------------------------------------------
# args: folds: a list of strings of the bracket notation folds for the seqs in seqsoral
#       seqsoral: a list of sequences (with the same list indeces as their structures in
#            the folds list
#       pos: list of integers; starting position list (like the inner lists from the
#            output of make_start_list)
#       le: integer; imposed length of the miR
# returns: 'left' or 'right'; a string indicating the side of the hairpin that the
#               candidate miR is on
# ------------------------------------------------------------------------------
# picks which side of the hairpin the miR is on based on the calculated secondary structures.
def pick_side(folds,seqsoral,pos,le):
    def pick_side_help(fold,seq,pos,le,r):
        def count_along(zp,zq,limit,step,seq):
            while zq<limit and zp<len(seq):
                if zp >=len(seq): print len(seq), zp
                if seq[zp]!='-': zq+=1
                zp+=step
            if  zp==len(seq): zp = zp - 1
            return zp
        def edge_helper(pos,p,el,er,r,le):
            a = pos - el - seq[el:pos].count('-') - 2
            b = er - p - seq[p:er].count('-')
            if a>b: return 'right'
            else: return 'left'
        brackets = { '(':1, '.':0, ')':-1, 'X':0 }
        p = count_along(pos,0,le,1,seq)
        way = sum(map(lambda a: brackets[a], fold[pos:p]))
        if way>0: return 'left'
        elif way<0: return 'right'
        else:
            sides = [filter(lambda a: a!='.', 'X'+fold[:pos])[-1],filter(lambda a: a!='.', fold[p:]+'X')[0]]
            if sides==['(',')']:
                el = pos-1
                while fold[el]=='.': el-=1
                er = p
                while fold[er]=='.': er+=1
                return edge_helper(pos,p,el,er,r,le)
            elif sides==['(','('] or sides==['X','(']: return 'left'
            elif sides==[')',')'] or sides==[')','X']: return 'right'
            elif sides==[')','(']:
                if fold[:pos].count('(')+fold[:pos].count(')') < fold[p:].count('(')+fold[p:].count(')'): return 'left'
                else: return 'right'
            elif sides==['X','X']: return edge_helper(0,len(fold),0,len(fold),r,le)
    answers = map(lambda r: pick_side_help(folds[r],seqsoral[r],pos[r],le,r), range(len(seqsoral)))
    if answers.count('left') >= answers.count('right'): return 'left'
    else: return 'right'



# mscore
# ------------------------------------------------------------------------------
# args: args: dictionary of arguements for funcitons
#       md: .matrix dictionary (output from parse_matrix)
# returns: out: dictionary of scores for each feature; keys refer to features (they
#               are the keys from fdict) and values are the assigned scores; the
#               additional key 'totscore' is bound to the sum of scores, and 'loc'
#               is bound to the integer of the position number in the first sequence
#               assuming that the first position is labelled "1".
# ------------------------------------------------------------------------------
# determines the mirscan score for a 21-mer at a particular position in the sequence
# referenced by ref#, position as indicated by the argument, by obtaining scores for
# each of the features in the matrix data structure return: dictionary with keys as
# names of score elements and values as the scores (floats)
def mscore(args,md,fdict,training=False):
    def totscore(lst):
        if training: return 'bogus_value'
        else: return sum(lst)
    def si(m,v):
        if training: return v
        else: return m[v]
    out = dict()
    for k in fdict.keys(): out[k] = si(md[k],fdict[k].ex(fdict[k],args))
    out['totscore'] = totscore(out.values())
    for n in range(len(args['orgs'])):
        out['loc_'+args['orgs'][n]] = args['pos'][n]+1
    return out




##############################################################
####  GENERATING ALIGNMENTS  #################################
##############################################################


def make_init_array(a,b):
    full = []
    for i in b:
        inner = []
        for j in a:
            inner.append(0)
        full.append(inner)
    return full

def make_init_scores(a,b,fgap,ngap=False):
    new = make_init_array(a,b)
    if ngap: new[0][0] = ngap-fgap
    for i in range(1,len(new)):
        new[i][0] = new[i-1][0]+fgap
    for j in range(1,len(new[0])):
        new[0][j] = new[0][j-1]+fgap
    return new

def make_init_trace(a,b):
    new = make_init_array(a,b)
    for i in range(1,len(new)):
        new[i][0] = 1
    for j in range(1,len(new[0])):
        new[0][j] = 2
    return new



# ------------------------------------------------------------------------------
# generates a local alignment of the two sequences (smith-waterman)
# ------------------------------------------------------------------------------
# seqa,seqb: sequences
# ngap,egap: open gap (new gap) and extend gap scores, respectively
# match,mis: match and mismatch scores, respectively
# fgap: score for having a gap at the edge of an alignment
# traceback: if given as True, returns a list with the score and traceback;
#            otherwise, the function just returns the score.
# NOTE: all scores should be the numbers that will be added to the
#       tally if the described feature is observed.  generally, 'match'
#       will be greater than zero and 'ngap', 'egap', 'mis', and 'fgap'
#       will be less than zero.
# ------------------------------------------------------------------------------
# if traceback==False:
#      returns the best local alignment score
# if traceback==True:
#      returns a 3-item list:
#      item 1: the the best local alignment score
#      item 2: the locally aligned portion of seqa, with gaps inserted
#      item 3: the locally aligned portion of seqb, with gaps inserted
# ------------------------------------------------------------------------------
def local_align(seqa,seqb,ngap,egap,match,mis,fgap,traceback=False):
    seqa = 'X'+seqa
    seqb = 'X'+seqb
    back = { "a":1, "b":2, "c":0 }
    ans = make_init_scores(seqa,seqb,0)
    aa = make_init_trace(seqa,seqb)
    max_score_coord = [0,0,0]
    for b in range(1,len(ans)):
        for a in range(1,len(ans[b])):
            c = dict()
            if a == len(ans[b])-1:
                c['a'] = max([ans[b-1][a]+fgap,0])
            elif aa[b-1][a] == 1:
                c['a'] = max([ans[b-1][a]+egap,0])
            else:
                c['a'] = max([ans[b-1][a]+ngap,0])
            if b == len(ans)-1:
                c['b'] = max([ans[b][a-1]+fgap,0])
            elif aa[b][a-1] == 2:
                c['b'] = max([ans[b][a-1]+egap,0])
            else:
                c['b'] = max([ans[b][a-1]+ngap,0])
            if seqa[a] == seqb[b]:
                c['c'] = max([ans[b-1][a-1]+match,0])
            else:
                c['c'] = max([ans[b-1][a-1]+mis,0])
            for i in 'cba':
                if c[i] == max(c.values()):
                    aa[b][a] = back[i]
            ans[b][a] = max(c.values())
            if ans[b][a] >= max_score_coord[0]:
                max_score_coord = [ans[b][a],b,a]
    if traceback:
        s = dict()
        n = dict()
        s[1] = [seqb,'',max_score_coord[1]]
        s[2] = [seqa,'',max_score_coord[2]]
        while ans[s[1][2]][s[2][2]]>0:
            if aa[s[1][2]][s[2][2]] == 0:
                for x in [1,2]:
                    s[x][1] = s[x][0][s[x][2]]+s[x][1]
                    s[x][2] += -1
            else:
                w = aa[s[1][2]][s[2][2]]
                s[w][1] = s[w][0][s[w][2]]+s[w][1]
                s[3-w][1] = '-'+s[3-w][1]
                s[w][2] += -1
        return [max_score_coord[0],s[2][1],s[1][1]]
    else: return max_score_coord[0]



# ------------------------------------------------------------------------------
# generates a global alignment of the two sequences (needleman-wunsch)
# ------------------------------------------------------------------------------
# seqa,seqb: sequences
# ngap,egap: open gap (new gap) and extend gap scores, respectively
# match,mis: match and mismatch scores, respectively
# fgap: score for having a gap at the edge of an alignment
# traceback: if given as True, returns a list with the score and traceback;
#            otherwise, the function just returns the score.
# NOTE: all scores should be the numbers that will be added to the
#       tally if the described feature is observed.  generally, 'match'
#       will be greater than zero and 'ngap', 'egap', 'mis', and 'fgap'
#       will be less than zero.
# ------------------------------------------------------------------------------
# if traceback==False:
#      returns the global alignment score
# if traceback==True:
#      returns a 3-item list:
#      item 1: the the global alignment score
#      item 2: the globally-aligned seqa, with gaps inserted
#      item 3: the globally-aligned seqb, with gaps inserted
# ------------------------------------------------------------------------------
def global_align(seqa,seqb,ngap,egap,match,mis,fgap,traceback=False):
    seqa = 'X'+seqa
    seqb = 'X'+seqb
    back = { "a":1, "b":2, "c":0 }
    if fgap: ans = make_init_scores(seqa,seqb,fgap)
    else: ans = make_init_scores(seqa,seqb,egap,ngap)
    aa = make_init_trace(seqa,seqb)
    for b in range(1,len(ans)):
        for a in range(1,len(ans[b])):
            c = dict()
            if fgap and a == len(ans[b])-1:
                c['a'] = ans[b-1][a]+fgap
            elif aa[b-1][a] == 1:
                c['a'] = ans[b-1][a]+egap
            else:
                c['a'] = ans[b-1][a]+ngap
            if fgap and b == len(ans)-1:
                c['b'] = ans[b][a-1]+fgap
            elif aa[b][a-1] == 2:
                c['b'] = ans[b][a-1]+egap
            else:
                c['b'] = ans[b][a-1]+ngap
            if seqa[a] == seqb[b]:
                c['c'] = ans[b-1][a-1]+match
            else:
                c['c'] = ans[b-1][a-1]+mis
            for i in 'cba':
                if c[i] == max(c.values()):
                    aa[b][a] = back[i]
            ans[b][a] = max(c.values())
    if traceback:
        s = dict()
        n = dict()
        s[1] = [seqb,'',len(seqb)-1]
        s[2] = [seqa,'',len(seqa)-1]
        while sum(map(lambda a: a[2], s.values())) > 0:
            if aa[s[1][2]][s[2][2]] == 0:
                for x in [1,2]:
                    s[x][1] = s[x][0][s[x][2]]+s[x][1]
                    s[x][2] += -1
            else:
                w = aa[s[1][2]][s[2][2]]
                s[w][1] = s[w][0][s[w][2]]+s[w][1]
                s[3-w][1] = '-'+s[3-w][1]
                s[w][2] += -1
        return [ans[len(seqb)-1][len(seqa)-1],s[2][1],s[1][1]]
    else: return ans[len(seqb)-1][len(seqa)-1]


# ------------------------------------------------------------------------------
# al: a 3-item list, where items 1 and 2 are the aligned sequences with gaps
# requires: items 1 and 2 are strings of the same length.
# returns: a string with a '*' or a ' ' (space) at each position in the alignment;
#          '*' if the identity is conserved, ' ' otherwise.
# ------------------------------------------------------------------------------
def make_star_line(al):
    newline = ''
    for n in range(len(al[1])):
        if al[1][n]==al[2][n]: newline+='*'
        else: newline+=' '
    return newline
