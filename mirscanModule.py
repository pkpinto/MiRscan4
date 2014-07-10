# coding: utf-8

import string, math, os, copy, sys

class Candidate:
    """A candidate is a potential miRNA hairpin or set of putatively orthologous
    miRNA hairpins.

    Attributes:
        name: Name of the candidate. It should be unique enough for the user to
            be able to distinguish between candidates based on the name alone.
    """
    def __init__(self, name, org_seq_dict):
        """
        Initialises the object with a name and dictionary of sequences.

        Args:
            name: A string unique enough for the user to be able to distinguish
                between candidates based on the name alone.
            org_seq_dict: a dictionary whose keys are organism keys (strings)
                and values are the corresponding miRNA hairpin sequences
                (strings of DNA letters; permitted characters are 'ACGTN').
        """
        self.name = name
        self._org_seq_dict = dict()
        for org in org_seq_dict.keys():
            self._org_seq_dict[org] = org_seq_dict[org]

    def organisms(self):
        """
        List of names (or a symbol representing the name/genome assembly, e.g.
        'dm2') of the organisms for which the candidate has candidate hairpins
        that are supposedly orthologous.
        """
        return self._org_seq_dict.keys()

    def seq(self, organism):
        """
        Given a organism name (organism), returns the genome sequences (DNA)
        corresponding to that candidate miRNA hairpin precursors.
        """
        return self._org_seq_dict[organism]


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
# from an instance of either the string_feature or number_feature class.  each of those
# classes has an internal state "type", which returns either 'string' or 'number' as
# a string, "kl" which is a list of possible key values from the scoring tables, and
# "kv", which is a method for taking a value returned by a "fx" method and turning it into
# an appropriate key value.  note that numerical values in kl must be in ascending order.
# there is also a function "pseudo" for each of the types of features which takes a string
# of counts and adds pseudocounts.
#
# string_feature
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
class string_feature:
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

# number_feature
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
class number_feature:
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

# File parsing  ################################################################

def parse_criteria(criteriafile):
    """
    Retrieves the mirscan criteria stored in a .py file.

    Any modules that the input file requires must be accessible from the
    working directory from which where this function is called.
    """
    with open(criteriafile) as input:
        code = input.read()
    criteria = dict()
    exec(code, criteria)
    return criteria

def get_background_files(trainfile):
    """
    Retrieves the background set filenames that are referred to in a .train
    file as a list of strings. If trainfile includes a relative path, it will
    added to the filename of the background file.
    """
    trainfile_path = os.path.dirname(trainfile)
    background_files = list()
    with open(trainfile) as input:
        for line in input:
            line = line.strip()
            if len(line) != 0 and line[0] == 'b':
                background_files.append(os.path.join(trainfile_path, line.split()[1]))
    return background_files

def parse_query(queryfile):
    """
    Detects the type of the query file and calls the correct parser.
    """
    if queryfile.split('.')[-1] == 'train':
        return parse_train(queryfile)
    elif queryfile.split('.')[-1] == 'fax':
        return parse_fax(queryfile)
    else:
        raise ValueError('Query file must be in \'.train\' or \'.fax\' format.')

def parse_train(trainfile, starts=False):
    """
    Retrieves the miRNA hairpin Candidates documented in a .train file.

    If starts==True, then two items are returned: the first is the (always
    returned) list of Candidates; the second is a list of dictionaries whose
    keys are org names and values are integers indicating the start position of
    the mature miRNA's 5p end in the corresponding hairpin sequence (indexed
    from 0). The index of a candidate in the first returned list matches the
    index of the start value dictionary in the second list.
    """
    uToT = string.maketrans('U','T')
    candidates = []
    startList = []
    firstOne = True
    with open(trainfile) as input:
        for line in input:
            columns = line.strip().split()
            if len(columns)!=0 and not(columns[0][0]=='#' or columns[0]=='b'):
                if columns[0]=='cn':
                    if firstOne:
                        firstOne = False
                    else:
                        candidates.append(Candidate(newName,newOrgToSeq))
                        if starts:
                            newStarts = dict()
                            for org in newOrgToSeq.keys():
                                newStarts[org] = newOrgToSeq[org].find(newOrgToMature[org])
                            startList.append(newStarts)
                    newName = columns[1]
                    newOrgToSeq = dict()
                    newOrgToMature = dict()
                elif columns[0]=='cm':
                    newOrgToMature[columns[1]] = columns[2].upper().translate(uToT,'')
                elif columns[0]=='ch':
                    newOrgToSeq[columns[1]] = columns[2].upper().translate(uToT,'')
    if firstOne:
        firstOne = False
    else:
        candidates.append(Candidate(newName,newOrgToSeq))
        if starts:
            newStarts = dict()
            for org in newOrgToSeq.keys():
                newStarts[org] = newOrgToSeq[org].find(newOrgToMature[org])
            startList.append(newStarts)

    __candidate_check(candidates)
    return (candidates,startList) if starts else candidates

def parse_fax(faxfile):
    """
    Retrieves the miRNA hairpin Candidates documented in a .fax file.
    """
    uToT = string.maketrans('U','T')
    candidates = []
    firstOne = True
    with open(faxfile) as input:
        for line in input:
            line = line.strip()
            if len(line)!=0 and line[0]=='#':
                line = ''
            if len(line)!=0:
                if line[0]=='>':
                    if firstOne:
                        firstOne = False
                    else:
                        candidates.append(Candidate(newName,newOrgToSeq))
                    newName = line[1:]
                    newOrgToSeq = dict()
                else:
                    columns = line.split()
                    newOrgToSeq[columns[0]] = columns[1].upper().translate(uToT,'')
    if firstOne:
        firstOne = False
    else:
        candidates.append(Candidate(newName, newOrgToSeq))

    __candidate_check(candidates)
    return candidates

def write_fax(candidates, faxfile):
    """
    Writes a set of miRNA hairpin Candidates to a .fax formated file. If
    filename=='stdout', writes to stdout.
    """
    __candidate_check(candidates)

    with (sys.stdout if faxfile.lower() == 'stdout' else open(faxfile, 'w')) as output:
        for c in candidates:
            output.write('>' + c.name + '\n')
            for org in c.organisms():
                output.write(org + '\t' + c.seq(org) + '\n')

def parse_matrix(matrixfile):
    """
    Retrieves the scoring matrix stored in a .matrix file.

    Returns:
        The scoring matrix as a dictionary of dictionaries. The outer
        dictionary's keys are feature names (i.e. the keys from fdict). The
        values (inner dictionaries) have keys which are elements of the keys
        lists (possible returned values) for that feature; the values are
        floating point numbers which are the scores corresponding to those
        particular feature values.
    """
    with open(matrixfile) as input:
        key = ''
        matrix = dict()
        first = 0
        for i in input.read().split('\n'):
            if not(i=='' or i[0]=='#'):
                if i==i.lstrip():
                    key = i.lstrip().rstrip()
                    matrix[key] = dict()
                else:
                    new = filter(lambda a: a!='', i.split('\t'))
                    matrix[key][new[0]] = float(new[1])
    return matrix

def __candidate_check(candidates):
    """
    The candidates from a file should each have hairpins from the same set of
    organisms, and those orgs should be referred to with the same key strings.
    """
    if len(candidates) > 0:
        initOrgs = candidates[0].organisms()
        if len(candidates) > 1:
            for c in candidates[1:]:
                for org in initOrgs:
                    if org not in c.organisms():
                        raise ValueError("candidate "+c.name+" is missing org: "+org)
                if len(c.organisms())!=len(initOrgs):
                    raise ValueError("candidate "+c.name+" has wrong number of orgs: "+str(c.organisms()))

def parse_scores(scoresfile, keys=False):
    """
    Retrieves the score sheets stroed in a .scr file. If keys are given, only
    those score features will be retrieved.
    """
    scores = list
    with open(scoresfile) as input:
        for line in input:
            columns = line.strip().split()
            if len(columns) > 1 and len(columns[0]) > 0 and columns[0][0] != '#':
                dd = dict()
                for n in range(len(columns[1:])/2):
                    dd[columns[2*n + 1]] = float(columns[2*n + 2])
                if keys:
                    d = {'name': columns[0]}
                    for k in keys:
                        d[k] = dd[k]
                    scores.append(d)
                else:
                    scores.append(dd)
    return scores

# General utilities  ###########################################################

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




# reverse_string
# ------------------------------------------------------------------------------
# args: st: a string
# returns: the same string with all characters appearing in the reverse order
# ------------------------------------------------------------------------------
# reverses the order of characters in a string
def reverse_string(st):
    li = []
    for i in st: li.append(i)
    li.reverse()
    return ''.join(li)




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
        fi,fo = os.popen2('RNAfold')
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



# make_bp_dicts
# ------------------------------------------------------------------------------
# args: fold: string which is a bracket-notation RNA secondary structure
# returns: bpl: a dictionary whose keys are integer positions in the secondary structure,
#               and whose values are 1 if the character in that positon is "(", 0 otherwise
#          bpr: a dictionary whose keys are integer positions in the secondary structure,
#               and whose values are 1 if the character in that positon is ")", 0 otherwise
#          bpd: a dictionary in which each paired position is a key and the position number
#               to which it is paired is its value
# ------------------------------------------------------------------------------
# takes bracket notation, makes dictionaries described.  the values in bpr and bpl are supposed
# to represent probabilities of those positions being paired.  since mfe structures are the
# input, probabilities are either 1 or 0.
def make_bp_dicts(fold):
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
    return bpl,bpr,bpd



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
