# coding: utf-8

import string, os, sys

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
    scores = list()
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
