#!/usr/bin/env python
# coding: utf-8

import math, types
import mirscanCriteria as msc

"""
The dictionary of all features considered for this criteria. All elements must
be an instance of the Feature class (or a child of it, eg. StringFeature).
"""
features = dict()

"""
The dictionary of DataItem which initialise required items of the cand_score
dictionary used in function Criteria.score.
Items defined here should not depend between themselves for the execution
order is not guaranteed. However, the standard items are always accessible and
can be made use of.
"""
data_items = dict()

"""
Default mature miRNA strand length for this criteria.
"""
mir_length = 22

# ------------------------------------------------------------------------------
"""
Requires: cand_data['hairpins']
"""

def _complexity_helper(seq,level_max):
    """
    args: seq: a string (a sequence of letters)
          level_max: an integer representing the longest repeated word that
                     will be sought
    returns: a normalization of the longest repeated word's length to the length of the
                overall sequence.  for sequence which is actually random, this value
                remains approximately constant with length of the sequence
    """
    tree = dict()
    sl = len(seq)
    pos = c = 0
    steps = []
    words = []
    while pos<sl:
        step = _find_in_tree(tree,1,pos,seq,level_max,sl)
        words.append(seq[pos:pos+step])
        pos += step
        steps.append(step)
        for n in range(1,level_max+2):
            if pos - n >= 0:
                _find_in_tree(tree,1,pos-n,seq,n,sl,1)
        c+=1
    longest = max(steps)
    return math.log(float(4**longest)/(len(seq)-longest*2)**2,4)

def _find_in_tree(tree,level,pos,seq,level_max,sl,build=False):
    """
    args: tree: a series of nested dictionaries representing all of the sequences
               that have been seen prior in the sequence, or a branch from that tree
               keys are nucleotide letters; values are dictionaries with key:value pairs
               for letters that appear after that letter (or series of letters) somewhere
               in the sequence
          level: integer; how many steps into the nested dictionary have been taken
          pos: integer; position in seq that is currently being compared to the tree
          seq: string; sequence being examined
          level_max: integer; maximum number of steps into the tree that will be taken
               (this is the same as the max repeated word length)
          sl: integer; length of 'seq'
          build: optional; if not false, updates (mutates) 'tree' to include the whole
               sequence up until the new position
    returns: integer corresponding to the number of steps that have been taken into the
                  tree (ie the number of steps forward that 'pos' should be advanced in
                  'seq' to get to the new starting position
    """
    if pos==sl: return level - 1
    elif tree.has_key(seq[pos]):
        if level==level_max: return level
        else: return _find_in_tree(tree[seq[pos]],level+1,pos+1,seq,level_max,sl,build)
    else:
        if build:
            tree[seq[pos]] = dict()
            if level!=level_max: return _find_in_tree(tree[seq[pos]],level+1,pos+1,seq,level_max,sl,build)
        if level==level_max: return level-1
        else: return max([level-1,1])

def _set_complexity(self, cand_data):
    return [('compl', map(lambda s: max([_complexity_helper(s, len(s) / 5),
                                         _complexity_helper(s[::-1], len(s) / 5)]),
                          cand_data['hairpins']))]

data_items['compl'] = msc.DataItem()
data_items['compl'].fx = types.MethodType(_set_complexity, data_items['compl'])

"""
Requires: cand_data['compl']
"""

def _get_complexity(self, cand_data):
    return sum(cand_data['compl']) / len(cand_data['compl'])

features['compl'] = msc.NumericalFeature()
features['compl'].fx = types.MethodType(_get_complexity, features['compl'])
features['compl'].bins = range(-2, 6)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
Requires: cand_data['folds']
"""

def _make_bp_dicts(fold):
    """
    args: fold: string which is a bracket-notation RNA secondary structure
    returns: bpl: a dictionary whose keys are integer positions in the secondary structure,
                  and whose values are 1 if the character in that position is "(", 0 otherwise
             bpr: a dictionary whose keys are integer positions in the secondary structure,
                  and whose values are 1 if the character in that position is ")", 0 otherwise
             bpd: a dictionary in which each paired position is a key and the position number
                  to which it is paired is its value

    takes bracket notation, makes dictionaries described.  the values in bpr and bpl are supposed
    to represent probabilities of those positions being paired.  since mfe structures are the
    input, probabilities are either 1 or 0.
    """
    bpd = dict()
    bpl = dict()
    bpr = dict()
    pll = []
    for n in range(len(fold)):
        if fold[n]=='(':
            bpl[n] = 1
            bpr[n] = 0
            pll.append(n)
        elif fold[n]==')':
            bpd[n] = pll.pop()
            bpd[bpd[n]] = n
            bpr[n] = 1
            bpl[n] = 0
        else:
            bpr[n] = 0
            bpl[n] = 0
    return bpl, bpr, bpd

def _make_bulge_list(fold, bpd):
    """
    args: fold: string of the bracket-notation nucleic acid 2ndary structure
          bpd: base pair dictionary, as generated by _make_bp_dicts
    returns: bulges: a list of bulges as described below

    creates a list of bulges for the hairpin to be passed to function fx5 whenever it is
    called for that hairpin.  each element in the list is a bulge; each bulge is an instance
    of class bpoint (defined here), which defines the bulge as internal state "b" and the
    opposite side of the bulge as internal state "o".  internal states "b" and "o" are each
    instances of class b2 (defined here), which has internal state "bef" (position of the
    nucleotide before the opening of the bulge) and internal state "aft" (position of the
    nucleotide which is paired to close the bulge).  unpaired positions at the beginning and
    end of the sequence are not counted as bulges.  all bulges appear twice in the list (such
    that the internal states "b" and "o" are switched).  loops appear once, with their "bef"
    and "aft" values stored in the internal state "b".  internal state "o" carries the value
    of 0 as its "bef" and "aft" internal states in the case of loops.
    """
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
    for i in fold:
        bulges_help.append(bpoint())
    for i in bulgesl:
        new = b2()
        if min(i)>0:
            new.bef = min(i)-1
        if max(i)<foldlength:
            new.aft = max(i)
        if new.bef and new.aft:
            bulges_first.append(bpoint())
            bulges_first[len(bulges_first)-1].b = new
            for j in range(min(i),max(i)):
                bulges_help[j].b = new
    loop_opp = b2()
    loop_opp.bef = loop_opp.aft = 0
    for i in bulges_first:
        if i.b.bef!=bpd[i.b.aft]:
            if bpd[i.b.bef]>i.b.bef and bulges_help[bpd[i.b.bef]-1].b:
                i.o = bulges_help[bpd[i.b.bef]-1].b
            elif bpd[i.b.bef]<i.b.bef and bulges_help[bpd[i.b.bef]-1].b:
                i.o = bulges_help[bpd[i.b.bef]-1].b
            else:
                i.o = b2()
                i.o.bef = bpd[i.b.aft]
                i.o.aft = bpd[i.b.bef]
        else:
            i.o = loop_opp
        if i.o and i.b.aft-i.b.bef+i.o.aft-i.o.bef>2:
            bulges.append(i)
    return bulges

def _set_bp(self, cand_data):
    bpl, bpr, bpd, bulges = list(), list(), list(), list()
    for fold in cand_data['folds']:
        _bpl, _bpr, _bpd = _make_bp_dicts(fold)
        bpl.append(_bpl)
        bpr.append(_bpr)
        bpd.append(_bpd)
        bulges.append(_make_bulge_list(fold, _bpd))
    return [('bpll', bpl),
            ('bprl', bpr),
            ('bpdl', bpd),
            ('bulgesl', bulges)]

data_items['bp'] = msc.DataItem()
data_items['bp'].fx = types.MethodType(_set_bp, data_items['bp'])

"""
Requires: cand_data['hairpins'], cand_data['bulgesl'], cand_data['mature_pos'],
            cand_data['length'], cand_data['alignment']
"""

def _bulge_sym_D_helper(seq,blg,pos,le):
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

def _get_bulge_sym_D(self, cand_data):
    """
    args: ar: arguements dictionary 'args' defined in mirscan.mirscan
    help args: seq: a sequence (or one sequence from an alignment)
               blg: a list of bulges as generated by _make_bulge_list for the 'al' entry of interest
               pos: integer starting position of the miR in the 'al' entry of interest
               le: integer imposed miR length
    help returns: integer number of unsymmetrically-bulged bases in the miR
    returns: arithmatic mean of help returned values over all candidate ortholog seqs

    for each bulge, counts number of unsymmetrically-bulged bases and sums them over all bulges
    uses individual structures
    """
    return float(sum(map(lambda r: _bulge_sym_D_helper(cand_data['hairpins'][r], cand_data['bulgesl'][r], cand_data['mature_pos'][r], cand_data['length']),
                         range(len(cand_data['alignment']))))) / len(cand_data['alignment'])

features['bulge_sym_D'] = msc.NumericalFeature()
features['bulge_sym_D'].fx = types.MethodType(_get_bulge_sym_D, features['bulge_sym_D'])
features['bulge_sym_D'].bins = range(10)

"""
Requires: cand_data['hairpins'], cand_data['folds'], cand_data['mature_pos'],
            cand_data['bpdl'], cand_data['alignment'], cand_data['mature_side'],
            cand_data['alignment']
"""

def _loop_dis_G_helper(fold,seq,bpd,pos,le,side,give_edges=False):
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
        while p-x>=pos and fold[p-x]!=')':
            x+=1
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
            while fold[p+x]!='(':
                x+=1
            if side=='left':
                return p,bpd[p+x]-(x+3)
            else:
                return edge_helper(pos,p,0,p+x,seq,le,side)
        def right_helper(pos,p,fold,seq,le,side):
            x = 1
            while fold[pos-x]!=')': x+=1
            if side=='right':
                if x-3<0:
                    return bpd[pos-x]+(1-x+3),pos
                else:
                    return bpd[pos-x]-(x-3),pos
            else:
                return edge_helper(pos,p,pos-x+1,len(seq),seq,le,side)
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
        elif sides==['X','(']:
            beginning,end = left_helper(pos,p,fold,seq,le,side)
        elif sides==[')','X']:
            beginning,end = right_helper(pos,p,fold,seq,le,side)
        elif sides==[')','(']:
            if fold[:pos].count('(')+fold[:pos].count(')') < fold[p:].count('(')+fold[p:].count(')'):
                beginning,end = left_helper(pos,p,fold,seq,le,side)
            else:
                beginning,end = right_helper(pos,p,fold,seq,le,side)
        elif sides==['X','X']:
            beginning,end = edge_helper(0,len(fold),0,len(fold),seq,le,side)
    if give_edges:
        return [beginning,end]
    else:
        return end - beginning

def _get_loop_dis_G(self, cand_data):
    """
    args: cand_data: arguements dictionary 'cand_data' defined in mirscan.mirscan
          give_edges: controls the returned value (optional; see below)
    help args: fold: string of the bracket-notation nucleic acid 2ndary structure
               seq: a sequence (or one sequence from an alignment)
               bpd: base pair dictionary, as generated by _make_bp_dicts
               pos: integer starting position of the miR in the 'al' entry of interest
               le: integer imposed miR length
               side: 'left' if the miR is on the left (5') side of the hairpin;
                     'right' if the miR is on the right (3') side of the hairpin;
                      as decided by pick_side
               give_edges: controls the returned value (optional; see below)
    help returns: if give_edges is False (not given), returns the length of the loop in nucleotides
                  (loop defined as the stretch between the miR and miR*).  otherwise, returns a two-item
                  list with the start and end coordinates of the loop
    returns: arithmatic mean of help returned values (loop lengths) over all candidate ortholog seqs
             if give_edges is not False, then it returns a list of two-item lists, where the
             outer list items are folds from cand_data['folds'] and the inner list has the start and
             end coordinates for that fold's loop (loop defined as the stretch between the miR and miR*)

    measures the loop size (loop defined as the stretch between the miR and miR*)
    using cand_data['mature_side'] and individual structures
    an implementation for counting loop length which is more ornate
    using cand_data['mature_side'] and individual structures
    """
    if self.give_edges:
        return map(lambda r: _loop_dis_G_helper(cand_data['folds'][r],cand_data['hairpins'][r],cand_data['bpdl'][r],cand_data['mature_pos'][r],cand_data['length'],cand_data['mature_side'],True), range(len(cand_data['alignment'])))
    else:
        return float(sum(map(lambda r: _loop_dis_G_helper(cand_data['folds'][r],cand_data['hairpins'][r],cand_data['bpdl'][r],cand_data['mature_pos'][r],cand_data['length'],cand_data['mature_side']), range(len(cand_data['alignment'])))))/len(cand_data['alignment'])

features['loop_dis_G'] = msc.NumericalFeature()
features['loop_dis_G'].fx = types.MethodType(_get_loop_dis_G, features['loop_dis_G'])
features['loop_dis_G'].bins = range(36)
features['loop_dis_G'].give_edges = False

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
Requires: cand_data['hairpins'], cand_data['mature_pos']
"""

def _get_nucX(self, cand_data):
    """
    args: cand_data: arguments dictionary 'cand_data' defined in mirscan.mirscan
               uses cand_data['hairpins'] and cand_data['mature_pos']
    returns: a one-letter string

    returns the one-letter string of the nucleotide at the position number of the sequence
    string provided.  'N' appears in genome sequence files and is also included in the nucs
    list.  if another character appears, it is printed and returned (an exception will result).
    """
    nuc = cand_data['hairpins'][0][ cand_data['mature_pos'][0] + self.position ]
    if nuc not in self.bins:
        raise ValueError('Unknown nucleotide present in sequence: ' + nuc)
    return nuc

features['nuc1_s1'] = msc.StringFeature()
features['nuc1_s1'].fx = types.MethodType(_get_nucX, features['nuc1_s1'])
features['nuc1_s1'].bins = ['A','T','C','G','N']
features['nuc1_s1'].position = 0

features['nuc9_s1'] = msc.StringFeature()
features['nuc9_s1'].fx = types.MethodType(_get_nucX, features['nuc9_s1'])
features['nuc9_s1'].bins = ['A','T','C','G','N']
features['nuc9_s1'].position = 8

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
Requires: cand_data['folds'], cand_data['mature_pos'], cand_data['mature_side'],
            cand_data['alignment']
"""

def _bp_matrix_A_helper(fold, pos, side):
    fd = { '(': 1, ')': 1, '.': 0 }
    if pos < 0 or pos >= len(fold):
        return 0
    elif side == 'left':
        return fd[fold[pos]]
    elif side == 'right':
        return fd[fold[pos]]

def _get_bp_matrix_A(self, cand_data):
    """
    Returns 'paired' if position is paired in all structures, otherwise the
    function returns 'unpaired'.
    """
    if float(sum(map(lambda r: _bp_matrix_A_helper(cand_data['folds'][r], self.position + cand_data['mature_pos'][r], cand_data['mature_side']),\
                     range(len(cand_data['alignment']))))) / len(cand_data['alignment']) == 1:
        return 'paired'
    else:
        return 'unpaired'

for n in range(mir_length):
    feature_name = 'bp_matrix_A_n' + str(n + 1)
    features[feature_name] = msc.StringFeature()
    features[feature_name].fx = types.MethodType(_get_bp_matrix_A, features[feature_name])
    features[feature_name].bins = ['paired', 'unpaired']
    features[feature_name].position = n

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
"""
Requires: cand_data['folds'], cand_data['mature_pos'], cand_data['mature_side'],
            cand_data['length'], cand_data['alignment']
"""

def _bp_matrix_C_helper(fold, pos, self_pos, side, le):
    fd = { '(': 1, ')': -1, '.': 0 }
    if side == 'left':
        if pos + self_pos < 0 or pos + self_pos >= len(fold):
            return 0
        else:
            return fd[fold[pos + self_pos]]
    elif side == 'right':
        if pos + le - 2 - self_pos >= len(fold) or pos + le - 2 - self_pos < 0:
            return 0
        else:
            return -fd[fold[pos + le - 2 - self_pos]]

def _get_bp_matrix_C(self, cand_data):
    """
    Similar to _get_bp_matrix_A except that when the mature strand is in the
    right (3') side of the hairpin, the sequence orientation is reversed so
    that the pairing probability is a function of position within the
    miRNA-miRNA* duplex.

    Returns 'paired' if position is paired in all structures, otherwise the
    function returns 'unpaired'.
    """
    if float(sum(map(lambda r: _bp_matrix_C_helper(cand_data['folds'][r], cand_data['mature_pos'][r], self.position, cand_data['mature_side'], cand_data['length']),\
                     range(len(cand_data['alignment']))))) / len(cand_data['alignment']) == 1:
        return 'paired'
    else:
        return 'unpaired'

for n in [-8, -7, -6, -5, -4, -3, -2, -1, 20, 21]:
    if n >= 0:
        feature_name = 'bp_matrix_C_n' + str(n + 1)
        features[feature_name] = msc.StringFeature()
        features[feature_name].fx = types.MethodType(_get_bp_matrix_C, features[feature_name])
        features[feature_name].bins = ['paired', 'unpaired']
        features[feature_name].position = n
    else:
        entry_name = 'bp_matrix_C_n' + str(n)
        features[feature_name] = msc.StringFeature()
        features[feature_name].fx = types.MethodType(_get_bp_matrix_C, features[feature_name])
        features[feature_name].bins = ['paired', 'unpaired']
        features[feature_name].position = n

# ------------------------------------------------------------------------------

"""
The Criteria class which will be exported from this file.
"""
criteria = msc.Criteria(mir_length, features, data_items)

# ------------------------------------------------------------------------------



# # method_con: conservtion matrix
# # ------------------------------------------------------------------------------
# # args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# # returns: all_match: string 'con' if the nuc identity at the positin given by
# #               self.pos is the same in all sequences; otherwise, returns 'non'
# # ------------------------------------------------------------------------------
# def method_con(self,ar,more_args=False):
# seqs = ar['hairpins']
# if more_args:
# pos = map(lambda n: n+more_args['start'], ar['mature_pos'])
# all_match = 0
# for n in range(more_args['length']):
#     add_this = 1
#     for n in range(1,len(seqs)):
#         if seqs[0][pos[0]]!=seqs[n][pos[n]]: add_this = 0
#     all_match += add_this
# return all_match
# else:
# pos = map(lambda n: n+self.pos, ar['mature_pos'])
# all_match = 'con'
# for n in range(1,len(seqs)):
#     if seqs[0][pos[0]]!=seqs[n][pos[n]]: all_match = 'non'
# return all_match
#
# for n in [1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20,21]:
# fdict['con'+str(n+1)] = msc.StringFeature()
# fdict['con'+str(n+1)].bins = ['con','non']
# fdict['con'+str(n+1)].pos = n
# fdict['con'+str(n+1)].fx = method_con
#
#
# # loop_vs_mir
# # uses method9_B (below) to figure out the fraction of the loop that is conserved;
# # uses method_con with more_args to get the fraction of the miR that is conserved;
# # returns the difference between those two fractions
# #
# # method343
# # ------------------------------------------------------------------------------
# # args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# # returns: float; difference between fraction of miR conserved and fraction of loop
# #               conserved; positive if miR is more conserved
# # ------------------------------------------------------------------------------
# # compares conservation of miR with conservation of loop; uses method_con
# def method343(self,ar):
# more_args = {'length':22, 'start':0}
# mircon = float(method_con(self,ar,more_args))/more_args['length']
# loopcon = method9_B(self,ar)
# return mircon - loopcon
# #
# # method9_B
# # ------------------------------------------------------------------------------
# # args: ar: arguements dictionary 'args' defined in mirscan.mirscan
# # returns: float; # conserved nucs / # total nucs in loop
# # ------------------------------------------------------------------------------
# # looks at loop conservation; uses method6_G
# def method9_B(self,ar):
# loop_list = method6_G(self,ar,1)
# seq_list = map(lambda r: ar['hairpins'][r][loop_list[r][0]:loop_list[r][1]+1],range(len(loop_list)))
# if max(map(len,seq_list))==0: return 0
# else:
# al = get_alignment(seq_list)
# match=0
# for n in range(len(al[0])):
#     if al[0][n]==al[1][n]: match+=1
# return float(match)/max(map(len,al))
#
# fdict['loop_vs_mir'] = msc.NumericalFeature()
# fdict['loop_vs_mir'].bins = map(lambda a: .05*a, range(-20,21))
# fdict['loop_vs_mir'].fx = method343
#
