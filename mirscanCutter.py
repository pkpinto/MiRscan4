#!/usr/bin/env python
# coding: utf-8

import sys, math, argparse
import mirscanIO as msio

def filtered_score(score, keys):
    """
    Computes total score as the sum of all scores identified in keys.
    """
    total = 0.0
    for k in keys:
        if k != 'name':
            total += score[k]
    return total

def compute_stats(mirs,kl):
    """
    Computes the mean score, stdev of scores, and min score calculated over
    filtered total score.
    """
    a = map(lambda m: filtered_score(m, kl), mirs)
    mean = sum(a) / len(a)
    dev = 0
    for n in a:
        diff = n - mean
        dev += diff * diff
    stdev = math.sqrt(dev / len(a))
    return mean, stdev, min(a)

def filter_scores(scores, keys, threshold):
    """
    Goes through a list of scores, gets a total score for keys defined in keys
    and returns those scores with a total above the threshold value.
    """
    selected = list()
    for s in scores:
        if filtered_score(s, keys) > threshold:
            selected.append(s)
    return selected


parser = argparse.ArgumentParser(description='''MiRscan3 Cutter.
                This script compares the score sheet files (.scr) of both the
                foreground and  background set of miRNA candidates. It computes
                the mean, standard deviation and minimum score of the
                foreground set to establish a threshold value. When the query
                file with background is provided (.fam of .fax), it
                proceeds to filter those candidates with a score above the
                threshold computed above.''',
            epilog='''Paulo Pinto, IEB-WWU, based on:
                http://bartellab.wi.mit.edu/softwareDocs/MiRscan3/Introduction.html''')

parser.add_argument(dest='fore_scorefile',
                    help='the foreground score sheet file (.scr)')
parser.add_argument(dest='back_scorefile',
                    help='the background score sheet file (.scr)')

parser.add_argument('--query_in', dest='in_queryfile',
                    help='the existing query file (.fam or .fax)')
parser.add_argument('--query_out', dest='out_queryfile', default='stdout',
                    help='the filtered query file (.fax) (default: stdout)')

args = parser.parse_args()

if args.fore_scorefile.split('.')[-1] != 'scr':
    raise ValueError('Foreground score sheet file must be in \'.scr\' format.')
if args.back_scorefile.split('.')[-1] != 'scr':
    raise ValueError('Background score sheet file must be in \'.scr\' format.')

if args.in_queryfile:
    if args.in_queryfile.split('.')[-1] != 'fam' and args.in_queryfile.split('.')[-1] != 'fax':
        raise ValueError('Query file must be in \'.fam\' or \'.fax\' format.')
    if args.out_queryfile.lower() != 'stdout' and args.out_queryfile.split('.')[-1] != 'fax':
        raise ValueError('Output query file must be in \'.fax\' format.')


# List of selected feature scores. ('name' key is always selected.)
score_keys = ['totscore']

fs = msio.parse_scores(args.fore_scorefile, score_keys)
bs = msio.parse_scores(args.back_scorefile, score_keys)

# The score distribution for the foreground set is analyzed to select a score
# threshold ('cut'), and the candidates which have scores above the threshold
# ('bcut' are those candidates' score dictionaries).
fmean, fstdev, fmin = compute_stats(fs, score_keys)
cut = fmin - fstdev / 2.0
bcut = filter_scores(bs, score_keys, cut)

print('\t'.join(['foreground:', 'mean=' + str(fmean), 'stdev=' + str(fstdev), 'cut=' + str(cut)]))
print('\t'.join(['background:', 'candidates above minimum=' + str(len(bcut)), 'total candidates=' + str(len(bs))]))

# If a source query file was provided, the passing candidates will be printed out.
if args.in_queryfile:
    candidates = msio.parse_query(args.in_queryfile)
    names = {c['name'] for c in bcut}
    selected_candidates = filter(lambda c: c.name in names, candidates)

    if len(selected_candidates) != len(bcut):
        raise ValueError('Number of scores and query candidates does not match.')
    msio.write_fax(selected_candidates, args.out_queryfile)
