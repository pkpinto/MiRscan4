# coding: utf-8

import string, os, random

class Feature(object):
    """
    Base class defining a computable feature of miRNA candidates. Usually a
    derived class (see below) should be used instead.

    This class defines a function (fx) used for computing the value associated
    with the feature. However, values are stored in bins keyed by values defined
    in a list of possible values (keys)

    Attributes:
        bins: List of allowed values for values of the feature. If the computed
              values are continuous, the bin function makes the translation to
              a discrete value.
    """
    def __init__(self):
        bins = None
    def fx(self, *args):
        """
        Function used to compute the value of the feature.
        """
        return None
    def pick_bin(self, val):
        """
        Bins the value of the feature (computed by fx) into the correct binned
        value (one of the values in bins).
        """
        return val
    def compute_bin(self, *args):
        """
        Shortcut method which computes the ferature value with fx and proceeds
        bin it with pick_bin.
        """
        return self.pick_bin(self.fx(*args))
    def smooth(self, count):
        """
        Given a set of values (count) counting the frequency of each bin value,
        this functions distributes the counts with the intent of smoothing the
        overall shape of the distribution.

        The default implementation does not really smooth the distribution, it
        just increments every bin by one. Though it effectivelly smoothens the
        distribution, its main objective is to prevent zero count bins which
        could cause divide-by-zero errors in processing.
        """
        return {b:(count[b] + 1) for b in count}

class StringFeature(Feature):
    def __init__(self):
        super(StringFeature, self).__init__()

class NumericalFeature(Feature):
    """
    Feature described by a numerical value. If it is a continuous value, bin
    discretises it into one of the values of bins.

    Attributes:
        bins: List of allowed values for values of the feature. They must be in
        ascending order. Requires that len(bins) >= 3
    """
    def __init__(self):
        super(NumericalFeature, self).__init__()
        self.smooth_round = 0
    def pick_bin(self, val):
        """
        Bins the value of the feature (computed by fx) into the correct binned
        value (one of the values in bins).

        The computed value val is placed in the bin where the value is less or
        equal than val and the value in the next bin is greater than val.
        """
        bin_idx = 1
        while bin_idx < len(self.bins) and self.bins[bin_idx] <= val:
            bin_idx += 1
        return self.bins[bin_idx - 1]
    def smooth(self, count):
        """
        Given a dict of values (count) counting the frequency of each bin value,
        this functions distributes the counts with the intent of smoothing the
        overall shape of the distribution.
        """
        if self.smooth_round == 2:
            self.smooth_round = 0
            return super(NumericalFeature, self).smooth(count)
        else:
            count_smooth = {b: 0 for b in count}
            kernel = {-1: 0.125, 0: 0.75, 1: 0.125}
            for n in range(1, len(self.bins) - 1):
                for x in kernel:
                    count_smooth[self.bins[n + x]] += kernel[x] * count[self.bins[n]]
            # left-most count
            count_smooth[self.bins[0]] += (1 - kernel[1]) * count[self.bins[0]]
            count_smooth[self.bins[1]] += kernel[1] * count[self.bins[0]]
            # right-most count
            count_smooth[self.bins[-1]] += (1 - kernel[-1]) * count[self.bins[-1]]
            count_smooth[self.bins[-2]] += kernel[-1] * count[self.bins[-1]]
            self.smooth_round += 1
            return self.smooth(count_smooth)

class DataItem:
    def fx(self, *args):
        """
        Function used to compute the value of the feature.
        """
        return None

class Criteria(object):
    """
    Class used for, given a dictionary of Features and DataItems, scoring sets
    of candidates miRNA sequences.
    """
    def __init__(self, mir_length, features, data_items=None):
        """Initialises the object with a dictionary of Features and DataItems.

        Args:
            mir_length: The Default mature miRNA strand length for this
                criteria.
            features: The dictionary of all features considered for this
                criteria. All elements must be an instance of the Feature class
                (or a child of it, eg. StringFeature).
            data_items: The dictionary of DataItem which initialise required
                items of the cand_score dictionary used in function
                Criteria.score. Items defined here should not depend between
                themselves for the execution order is not guaranteed. However,
                the standard items are always accessible and can be made use of.
        """
        self.mir_length = mir_length
        self.features = features
        self.data_items = data_items

    def score(self, candidates, matrix=None):
        """
        Scores a set of Candidates. Returns a list of scores, one for each
        candidate. The candidate score is computed from a dictionary of
        features. Each feature score is then compared to the pre-computed matrix
        reference value.

        If a scoring matrix is not provided, the function falls into training
        mode. This allows for a scoring matrix to be created. Each feature score
        is simply returned for further processing (see mirscanTrainer).
        """
        cand_score = list()
        for cand in candidates:
            # Candidate data dictionary stores a variety of variables useful for
            # assessing the score of this candidate
            cand_data = dict()
            cand_data['length'] = self.mir_length
            # Ensure that the list of organisms is identical across candidates
            cand_data['organisms'] = cand.organisms()
            cand_data['organisms'].sort()

            # sequences
            cand_data['hairpins'] = map(lambda org: cand.hairpin(org), cand_data['organisms'])
            # folds
            if not cand.has_folds():
                cand.compute_folds()
            cand_data['folds'] = map(lambda org: cand.fold(org), cand_data['organisms'])
            # alignments
            if len(cand.organisms()) > 1:
                # cand_data['alignment'] = get_alignment(args['hairpins'])
                cand_data['alignment'] = cand_data['hairpins'][0]
            else:
                cand_data['alignment'] = cand_data['hairpins']

            # Extra entries to cand_data defined in the criteria file
            for ditem in self.data_items:
                for k,v in self.data_items[ditem].fx(cand_data):
                    cand_data[k] = v

            # Get the scores for all of the possible (or specified) start positions
            candlist = []
            if cand.has_matures():
                cand_data['mature_pos'] = map(lambda org: cand.hairpin(org).find(cand.mature(org)), cand_data['organisms'])
                cand_data['mature_side'] = self._pick_side(cand_data['folds'], cand_data['hairpins'], cand_data['mature_pos'], cand_data['length'])
                candlist.append(self._cand_score(cand_data, matrix))
            else:
                start_list = self._make_start_list(cand_data['hairpins'], cand_data['alignment'], cand_data['length'])
                # training mode: pick just one start point
                if not matrix:
                    start_list = [random.choice(start_list)]
                for i in start_list:
                    cand_data['mature_pos'] = i
                    cand_data['mature_side'] = self._pick_side(cand_data['folds'], cand_data['hairpins'], cand_data['mature_pos'], cand_data['length'])
                    candlist.append(self._cand_score(cand_data, matrix))

            # sort the start position options by the totscores and return the highest
            npl = [(x['total_score'],x) for x in candlist]
            npl.sort()
            candlist = [val for (key, val) in npl]
            bestCand = candlist[-1]
            bestCand['name'] = cand.name
            cand_score.append(bestCand)

        return cand_score

    def _cand_score(self, cand_data, matrix=None):
        """
        args: cand_data: dictionary of arguements for funcitons
              matrix: .matrix dictionary (output from parse_matrix)
        returns: score: dictionary of scores for each feature; keys refer to features (they
                      are the keys from fdict) and values are the assigned scores; the
                      additional key 'total_score' is bound to the sum of scores, and 'mature_pos'
                      is bound to the integer of the position number in the first sequence
                      assuming that the first position is labelled "1".

        determines the mirscan score for a 21-mer at a particular position in the sequence
        referenced by ref#, position as indicated by the argument, by obtaining scores for
        each of the features in the matrix data structure return: dictionary with keys as
        names of score elements and values as the scores (floats)
        """
        score = dict()
        for feature in self.features:
            if matrix:
                score[feature] = matrix[feature][ self.features[feature].compute_bin(cand_data) ]
            else:
                score[feature] = self.features[feature].fx(cand_data)
        score['total_score'] = sum(score.values()) if matrix else -1
        for i,org in enumerate(cand_data['organisms']):
            score['mature_pos_' + org] = cand_data['mature_pos'][i] + 1
        return score

    def _make_start_list(self,seq_list,al,length):
        """
        args: seq_list: a list of the query sequences with no gaps ('-' characters)
              al: a list of the query sequences as they would appear in an alignment
                   with one another (ie with gaps, all strings should be the same length)
              length: integer length of a miR
        returns: pos_list: a list of lists.  each element in the outer list corresponds
                      to a set of start positions for a miR.  the first n elements of the
                      the inner list gives the starting position in each of the n sequences
                      of seq_list, in the same order as the sequences they refer to.  the
                      n+1 element gives the starting position in the alignment.

        generates a list of lists, where the elements of the list are possible start positions for the mir
        in each of the sequences; start positions in the returned list are in the order of the sequences
        within the input sequence list.
        """
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

    def _pick_side(self,folds,seqsoral,pos,le):
        """
        args: folds: a list of strings of the bracket notation folds for the seqs in seqsoral
              seqsoral: a list of sequences (with the same list indeces as their structures in
                   the folds list
              pos: list of integers; starting position list (like the inner lists from the
                   output of _make_start_list)
              le: integer; imposed length of the miR
        returns: 'left' or 'right'; a string indicating the side of the hairpin that the
                      candidate miR is on

        picks which side of the hairpin the miR is on based on the calculated secondary structures.
        """
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



# ##############################################################
# ####  GENERATING ALIGNMENTS  #################################
# ##############################################################
#
# get_alignment
# ------------------------------------------------------------------------------
# args: two sequences
# returns: alignment as list of strings (input sequences + gaps)
# ------------------------------------------------------------------------------
# finds alignment of two sequences
# def get_alignment(sa):
#     return global_align(sa[0],sa[1],-8,-3,1,-3,-2,'get traceback')[1:]
#
#
# def make_init_array(a,b):
#     full = []
#     for i in b:
#         inner = []
#         for j in a:
#             inner.append(0)
#         full.append(inner)
#     return full
#
# def make_init_scores(a,b,fgap,ngap=False):
#     new = make_init_array(a,b)
#     if ngap: new[0][0] = ngap-fgap
#     for i in range(1,len(new)):
#         new[i][0] = new[i-1][0]+fgap
#     for j in range(1,len(new[0])):
#         new[0][j] = new[0][j-1]+fgap
#     return new
#
# def make_init_trace(a,b):
#     new = make_init_array(a,b)
#     for i in range(1,len(new)):
#         new[i][0] = 1
#     for j in range(1,len(new[0])):
#         new[0][j] = 2
#     return new
#
#
# # ------------------------------------------------------------------------------
# # generates a local alignment of the two sequences (smith-waterman)
# # ------------------------------------------------------------------------------
# # seqa,seqb: sequences
# # ngap,egap: open gap (new gap) and extend gap scores, respectively
# # match,mis: match and mismatch scores, respectively
# # fgap: score for having a gap at the edge of an alignment
# # traceback: if given as True, returns a list with the score and traceback;
# #            otherwise, the function just returns the score.
# # NOTE: all scores should be the numbers that will be added to the
# #       tally if the described feature is observed.  generally, 'match'
# #       will be greater than zero and 'ngap', 'egap', 'mis', and 'fgap'
# #       will be less than zero.
# # ------------------------------------------------------------------------------
# # if traceback==False:
# #      returns the best local alignment score
# # if traceback==True:
# #      returns a 3-item list:
# #      item 1: the the best local alignment score
# #      item 2: the locally aligned portion of seqa, with gaps inserted
# #      item 3: the locally aligned portion of seqb, with gaps inserted
# # ------------------------------------------------------------------------------
# def local_align(seqa,seqb,ngap,egap,match,mis,fgap,traceback=False):
#     seqa = 'X'+seqa
#     seqb = 'X'+seqb
#     back = { "a":1, "b":2, "c":0 }
#     ans = make_init_scores(seqa,seqb,0)
#     aa = make_init_trace(seqa,seqb)
#     max_score_coord = [0,0,0]
#     for b in range(1,len(ans)):
#         for a in range(1,len(ans[b])):
#             c = dict()
#             if a == len(ans[b])-1:
#                 c['a'] = max([ans[b-1][a]+fgap,0])
#             elif aa[b-1][a] == 1:
#                 c['a'] = max([ans[b-1][a]+egap,0])
#             else:
#                 c['a'] = max([ans[b-1][a]+ngap,0])
#             if b == len(ans)-1:
#                 c['b'] = max([ans[b][a-1]+fgap,0])
#             elif aa[b][a-1] == 2:
#                 c['b'] = max([ans[b][a-1]+egap,0])
#             else:
#                 c['b'] = max([ans[b][a-1]+ngap,0])
#             if seqa[a] == seqb[b]:
#                 c['c'] = max([ans[b-1][a-1]+match,0])
#             else:
#                 c['c'] = max([ans[b-1][a-1]+mis,0])
#             for i in 'cba':
#                 if c[i] == max(c.values()):
#                     aa[b][a] = back[i]
#             ans[b][a] = max(c.values())
#             if ans[b][a] >= max_score_coord[0]:
#                 max_score_coord = [ans[b][a],b,a]
#     if traceback:
#         s = dict()
#         n = dict()
#         s[1] = [seqb,'',max_score_coord[1]]
#         s[2] = [seqa,'',max_score_coord[2]]
#         while ans[s[1][2]][s[2][2]]>0:
#             if aa[s[1][2]][s[2][2]] == 0:
#                 for x in [1,2]:
#                     s[x][1] = s[x][0][s[x][2]]+s[x][1]
#                     s[x][2] += -1
#             else:
#                 w = aa[s[1][2]][s[2][2]]
#                 s[w][1] = s[w][0][s[w][2]]+s[w][1]
#                 s[3-w][1] = '-'+s[3-w][1]
#                 s[w][2] += -1
#         return [max_score_coord[0],s[2][1],s[1][1]]
#     else: return max_score_coord[0]
#
#
# # ------------------------------------------------------------------------------
# # generates a global alignment of the two sequences (needleman-wunsch)
# # ------------------------------------------------------------------------------
# # seqa,seqb: sequences
# # ngap,egap: open gap (new gap) and extend gap scores, respectively
# # match,mis: match and mismatch scores, respectively
# # fgap: score for having a gap at the edge of an alignment
# # traceback: if given as True, returns a list with the score and traceback;
# #            otherwise, the function just returns the score.
# # NOTE: all scores should be the numbers that will be added to the
# #       tally if the described feature is observed.  generally, 'match'
# #       will be greater than zero and 'ngap', 'egap', 'mis', and 'fgap'
# #       will be less than zero.
# # ------------------------------------------------------------------------------
# # if traceback==False:
# #      returns the global alignment score
# # if traceback==True:
# #      returns a 3-item list:
# #      item 1: the the global alignment score
# #      item 2: the globally-aligned seqa, with gaps inserted
# #      item 3: the globally-aligned seqb, with gaps inserted
# # ------------------------------------------------------------------------------
# def global_align(seqa,seqb,ngap,egap,match,mis,fgap,traceback=False):
#     seqa = 'X'+seqa
#     seqb = 'X'+seqb
#     back = { "a":1, "b":2, "c":0 }
#     if fgap: ans = make_init_scores(seqa,seqb,fgap)
#     else: ans = make_init_scores(seqa,seqb,egap,ngap)
#     aa = make_init_trace(seqa,seqb)
#     for b in range(1,len(ans)):
#         for a in range(1,len(ans[b])):
#             c = dict()
#             if fgap and a == len(ans[b])-1:
#                 c['a'] = ans[b-1][a]+fgap
#             elif aa[b-1][a] == 1:
#                 c['a'] = ans[b-1][a]+egap
#             else:
#                 c['a'] = ans[b-1][a]+ngap
#             if fgap and b == len(ans)-1:
#                 c['b'] = ans[b][a-1]+fgap
#             elif aa[b][a-1] == 2:
#                 c['b'] = ans[b][a-1]+egap
#             else:
#                 c['b'] = ans[b][a-1]+ngap
#             if seqa[a] == seqb[b]:
#                 c['c'] = ans[b-1][a-1]+match
#             else:
#                 c['c'] = ans[b-1][a-1]+mis
#             for i in 'cba':
#                 if c[i] == max(c.values()):
#                     aa[b][a] = back[i]
#             ans[b][a] = max(c.values())
#     if traceback:
#         s = dict()
#         n = dict()
#         s[1] = [seqb,'',len(seqb)-1]
#         s[2] = [seqa,'',len(seqa)-1]
#         while sum(map(lambda a: a[2], s.values())) > 0:
#             if aa[s[1][2]][s[2][2]] == 0:
#                 for x in [1,2]:
#                     s[x][1] = s[x][0][s[x][2]]+s[x][1]
#                     s[x][2] += -1
#             else:
#                 w = aa[s[1][2]][s[2][2]]
#                 s[w][1] = s[w][0][s[w][2]]+s[w][1]
#                 s[3-w][1] = '-'+s[3-w][1]
#                 s[w][2] += -1
#         return [ans[len(seqb)-1][len(seqa)-1],s[2][1],s[1][1]]
#     else: return ans[len(seqb)-1][len(seqa)-1]
#
#
# # ------------------------------------------------------------------------------
# # al: a 3-item list, where items 1 and 2 are the aligned sequences with gaps
# # requires: items 1 and 2 are strings of the same length.
# # returns: a string with a '*' or a ' ' (space) at each position in the alignment;
# #          '*' if the identity is conserved, ' ' otherwise.
# # ------------------------------------------------------------------------------
# def make_star_line(al):
#     newline = ''
#     for n in range(len(al[1])):
#         if al[1][n]==al[2][n]: newline+='*'
#         else: newline+=' '
#     return newline
