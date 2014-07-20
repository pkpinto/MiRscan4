# coding: utf-8

import string, os, sys

class Candidate(object):
    """A candidate is a potential miRNA hairpin or set of putatively orthologous
    miRNA hairpins.

    Attributes:
        name: Name of the candidate. It should be unique enough for the user to
            be able to distinguish between candidates based on the name alone.
    """
    def __init__(self, name, org_hairpin_dict, org_mature_dict=None, org_fold_dict=None):
        """
        Initialises the object with a name and dictionary of sequences.

        Args:
            name: A string unique enough for the user to be able to distinguish
                between candidates based on the name alone.
            org_hairpin_dict: a dictionary whose keys are organism keys (strings)
                and values are the corresponding DNA sequence of the miRNA
                hairpin sequences (strings of DNA letters; permitted characters
                are 'ATCGN').
            org_mature_dict: a dictionary whose keys are organism keys (strings)
                and values are the corresponding DNA sequence of the miRNA
                mature strand sequences (strings of RNA letters; permitted
                characters are 'ATCGN'). Mature strand sequences to all
                organisms present in org_hairpin_dict must be provided,
                additional organisms not present in org_hairpin_dict are ignored.
            org_fold_dict: a dictionary whose keys are organism keys (strings)
                and values are the fold of the miRNAs in dot-bracket notation.
                Folds for all organisms present in org_hairpin_dict must be
                provided, additional organisms not present in org_hairpin_dict
                are ignored.
            """
        self.name = name
        self._org_hairpin_dict = dict()
        for org in org_hairpin_dict:
            self._org_hairpin_dict[org] = org_hairpin_dict[org]
        if org_mature_dict:
            self._org_mature_dict = dict()
            # ensure all organisms for which there is a harpin also have a
            # mature strand assigned, but only those (extra mature sequences
            # are ignored)
            for org in self._org_hairpin_dict:
                self._org_mature_dict[org] = org_mature_dict[org]
        else:
            self._org_mature_dict = None
        if org_fold_dict:
            self._org_fold_dict = dict()
            for org in self._org_fold_dict:
                self._org_fold_dict[org] = org_fold_dict[org]
        else:
            self._org_fold_dict = None

    def organisms(self):
        """
        List of names (or a symbol representing the name/genome assembly, e.g.
        'dm2') of the organisms for which the candidate has candidate hairpins
        that are supposedly orthologous.
        """
        return self._org_hairpin_dict.keys()

    def hairpin(self, organism):
        """
        Given a organism name (organism), returns the DNA sequence of the
        candidate miRNA hairpin precursor.
        """
        return self._org_hairpin_dict[organism]

    def has_matures(self):
        """
        Returns True when mature sequences are available
        """
        return True if self._org_mature_dict else False
    def mature(self, organism):
        """
        Given a organism name (organism), returns the DNA sequence of the
        candidate miRNA mature strand. Returns None if no mature strand
        sequences were provided on initialisation.
        """
        if self._org_mature_dict:
            return self._org_hairpin_dict[organism]
        else:
            return None

    def fold(self, organism):
        """
        Given a organism name (organism), returns the corresponding fold in
        dot-bracket notation. Returns None if the folds are not available.
        """
        if self._org_fold_dict:
            return self._org_fold_dict[organism]
        else:
            return None
    def has_folds(self):
        """
        Returns True when fold information is available
        """
        return True if self._org_fold_dict else False

    def compute_folds(self):
        """
        RNAfold (from the ViennaRNA package) must be executable from the command
        line by typing 'RNAfold'.
        """
        fold_cmd = 'RNAfold --noPS'
        self._org_fold_dict = dict()
        for org in self._org_hairpin_dict:
            fi,fo = os.popen2(fold_cmd)
            fi.write(self._org_hairpin_dict[org])
            fi.close()
            fold = fo.read().split('\n')[1].split(' ')[0]
            self._org_fold_dict[org] = fold

def parse_criteria(criteriafile):
    """
    Retrieves the mirscan criteria stored in a .py file.

    Any modules that the input file requires must be accessible from the
    working directory from which where this function is called.
    """
    with open(criteriafile, 'r') as input:
        code = input.read()
    environment = dict()
    exec(code, environment)
    return environment['criteria']

def parse_query(queryfile):
    """
    Detects the type of the query file and calls the correct parser.
    """
    if queryfile.split('.')[-1] == 'fam':
        return parse_fam(queryfile)
    elif queryfile.split('.')[-1] == 'fax':
        return parse_fax(queryfile)
    else:
        raise ValueError('Query file must be in \'.fam\' or \'.fax\' format.')

def parse_fam(famfile, starts=False):
    """
    Retrieves the miRNA hairpin and mature Candidates documented in a .fam file.
    """
    RNAtoDNA = string.maketrans('U', 'T')
    candidates = list()
    open_block = False
    with open(famfile, 'r') as input:
        for line in input:
            line = line.strip()
            if len(line) == 0:
                continue
            elif line[0] == '#':
                continue
            else:
                if line[0] == '>':
                    if open_block:
                        candidates.append(Candidate(__name, __org_hairpin_dict, __org_mature_dict))
                        # open_block = False
                    open_block = True
                    __name = line[1:].strip()
                    __org_hairpin_dict = dict()
                    __org_mature_dict = dict()
                else:
                    strand,org,seq = line.split()
                    if strand == 'h':
                        __org_hairpin_dict[org] = seq.upper().translate(RNAtoDNA, '')
                    elif strand == 'm':
                        __org_mature_dict[org] = seq.upper().translate(RNAtoDNA, '')
    if open_block:
        candidates.append(Candidate(__name, __org_hairpin_dict, __org_mature_dict))
        # open_block = False
    return candidates

def parse_fax(faxfile):
    """
    Retrieves the miRNA hairpin Candidates documented in a .fax file.
    """
    RNAtoDNA = string.maketrans('U', 'T')
    candidates = list()
    open_block = False
    with open(faxfile, 'r') as input:
        for line in input:
            line = line.strip()
            if len(line) == 0:
                continue
            elif line[0] == '#':
                continue
            else:
                if line[0] == '>':
                    if open_block:
                        candidates.append(Candidate(__name, __org_hairpin_dict))
                        # open_block = False
                    open_block = True
                    __name = line[1:].strip()
                    __org_hairpin_dict = dict()
                else:
                    org,seq = line.split()
                    __org_hairpin_dict[org] = seq.upper().translate(RNAtoDNA, '')
    if open_block:
        candidates.append(Candidate(__name, __org_hairpin_dict))
        # open_block = False
    __check_candidates(candidates)
    return candidates

def write_fax(candidates, faxfile):
    """
    Writes a set of miRNA hairpin Candidates to a .fax formated file. If
    faxfile=='stdout', writes to stdout.
    """
    __check_candidates(candidates)

    with (sys.stdout if faxfile.lower() == 'stdout' else open(faxfile, 'w')) as output:
        for c in candidates:
            output.write('>' + c.name + '\n')
            for org in c.organisms():
                output.write(org + ' ' + c.hairpin(org) + '\n')

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
    with open(matrixfile, 'r') as input:
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

def __check_candidates(candidates):
    """
    The candidates from a file should each have hairpins from the same set of
    organisms, and those orgs should be referred to with the same key strings.
    """
    if len(candidates) > 0:
        organism_list = candidates[0].organisms()
        if len(candidates) > 1:
            for c in candidates[1:]:
                for org in organism_list:
                    if org not in c.organisms():
                        raise ValueError('Candidate %s should contain a sequence for organism %s.' % (c.name, org))
                if len(c.organisms()) != len(organism_list):
                    raise ValueError('Candidate %s does not have sequences for the expected number of organisms.' % c.name)

def parse_scores(scoresfile, keys=False):
    """
    Retrieves the score sheets stored in a .scr file. If keys are given, only
    those score features will be retrieved.
    """
    with open(scoresfile, 'r') as input:
        unselected_scores = eval(input.read())
    if keys:
        scores = list()
        for candidate_score in unselected_scores:
            s = {'name': candidate_score['name']}
            for k in keys:
                s[k] = candidate_score[k]
            scores.append(s)
        return scores
    else:
        return unselected_scores
