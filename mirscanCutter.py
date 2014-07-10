#!/usr/bin/env python
# coding: utf-8

# This script collects the scores calculated for the foreground cadidates and
# finds the mean and standard deviation of the foreground score distribution.
# It decides where to cut the background: the point where the score which is
# one standard deviation below the lowest score from the foreground set.

import sys, math, argparse
import mirscanIO as msio

def filtered_score(score, keys):
    """
    Computes total score as the sum of all scores identified in keys.
    """
    total = 0
    for k in keys:
        total += score[k]
    return total

def compute_stats(mirs,kl):
    """
    Computes the mean score, stdev of scores, and min score calculated over
    filtered total score.
    """
    a = map(lambda m: filtered_score(m,kl), mirs)
    mean = sum(a)/len(a)
    dev = 0
    for n in a:
        diff = n-mean
        dev += diff*diff
    stdev = math.sqrt(dev/len(a))
    return mean,stdev,min(a)

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


parser = argparse.ArgumentParser(description='MiRscan3 Cutter',
            epilog='Paulo Pinto, IEB-WWU, based on:\nhttp://bartellab.wi.mit.edu/softwareDocs/MiRscan3/Introduction.html')

parser.add_argument(dest='fore_scorefile',
                    help='the foreground score sheet file (.scr)')
parser.add_argument(dest='back_scorefile',
                    help='the background score sheet file (.scr)')

parser.add_argument('--query_in', dest='in_queryfile',
                    help='the existing query file (.train or .fax)')
parser.add_argument('--query_out', dest='out_queryfile', default='stdout',
                    help='the filtered query file (.fax) (default: stdout)')

args = parser.parse_args()

if args.fore_scorefile.split('.')[-1] != 'scr':
    raise ValueError('Foreground score sheet file must be in \'.scr\' format.')
if args.back_scorefile.split('.')[-1] != 'scr':
    raise ValueError('Background score sheet file must be in \'.scr\' format.')

if args.in_queryfile:
    if args.in_queryfile.split('.')[-1] != 'train' and args.in_queryfile.split('.')[-1] != 'fax':
        raise ValueError('Query file must be in \'.train\' or \'.fax\' format.')
    if args.out_queryfile.lower() != 'stdout' and args.out_queryfile.split('.')[-1] != 'fax':
        raise ValueError('Output query file must be in \'.fax\' format.')


# this list can be modified to change the behavior of score_cut.py, but here,
# it is set to use the pre-computed (by mirscan) sum of all of the individual
# feature scores, referred to as 'totscore' in mirscan's output.
score_keys = ['totscore']

fs = msio.parse_scores(args.fore_scorefile, score_keys)
bs = msio.parse_scores(args.back_scorefile, score_keys)

# the score distribution for the foreground set is analyzed to select a score
# threshold ('cut'), and the candidates are filtered for having scores above
# the threshold ('bcut' are those candidates' score dictionaries).
fmean,fstdev,fmin = compute_stats(fs, score_keys)
cut = fmin - fstdev / 2.0
bcut = filter_scores(bs, score_keys, cut)

print('\t'.join(['foreground:', 'mean=' + str(fmean), 'stdev=' + str(fstdev), 'cut=' + str(cut)]))
print('\t'.join(['background:', 'candidates above minimum=' + str(len(bcut)), 'total candidates=' + str(len(bs))]))

# if the appropriate arguments were provided, the passing candidates will be
# placed into a new .fax file.
if args.in_queryfile:
    candidates = msio.parse_query(args.in_queryfile)
    names = {c['name'] for c in bcut}
    selected_candidates = filter(lambda c: c.name in names, candidates)

    if len(selected_candidates) != len(bcut):
        raise ValueError('number of scored (' + str(len(bcut)) +') and queried (' \
                        + str(len(selected_candidates)) + ') candidates don\'t match after filtering.')
    msio.write_fax(selected_candidates, args.out_queryfile)
