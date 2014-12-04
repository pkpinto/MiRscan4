MiRscan4
========

Fork of MiRscan3, fixing bugs and updating the code.

MiRscan was originally developed for the detection of miRNAs in C. Elegans [Lim
et al. 2003a] and vertebrates [Lim et al. 2003b]. More recently an updated
version has been release and applied to the detection of Drosophila miRNAs
[Ruby et al. 2007].

### mirscanTrainer.py
This script takes a set of foreground (.fam, the training procedure requires
mature strand sequences) and background (.fam or .fax) miRNA candidates, a
mirscan criteria file (.py) with rules for their evaluation and outputs a
scoring matrix file (.matrix) which can be used for scoring additional miRNAs.

### mirscanScorer.py
Given a query file listing miRNA candidates (.fam or .fax), this script scores
them using an (also given) mirscan criteria file (.py) and scoring matrix file
(.matrix). The output is done to a score sheet file (.scr).

### mirscanCutter.py
This script compares the score sheet files (.scr) of both the foreground and
background set of miRNA candidates. It computes the mean, standard deviation
and minimum score of the foreground set to establish a threshold value.

When the query file with background is provided (.fam of .fax), it proceeds to
filter those candidates with a score above the threshold computed above.

## File formats

### .fax
Based on the fasta format, this file format can store several (putatively)
orthologous sequences grouping together the organisms where their are present.

```
>mir-2a-2
dm2 ATCTAAGCCTCATCAAGTGGTTGTGATATGGATACCCAACGCATATCACAGCCAGCTTTGATGAGCTAGGAT
dp3 AUCUAAGCCUCAUCAAGUGGUUGUGAUAUGGAUACCCAACGCAUAUCACAGCCAGCUUUGAUGAGCUAGGAU
>mir-2b-1
ch dm2 CTTCAACTGTCTTCAAAGTGGCAGTGACATGTTGTCAACAATATTCATATCACAGCCAGCTTTGAGGAGCGTTGCGG
ch dp3 CUGCGACGCUCUUUAAAGUGGCGGUGACGUGUUGGUAAUAAUAUUCAUAUCACAGCCAGCUUUGAGGAGCGUUGCGG
```

### .fam
Also based on fasta, .fam files can additionally store, for each species, a
hairpin and mature strand sequence.
```
>mir-2a-2
m dm2 uaucacagccagcuuugaugagc
h dm2 ATCTAAGCCTCATCAAGTGGTTGTGATATGGATACCCAACGCATATCACAGCCAGCTTTGATGAGCTAGGAT
m dp3 uaucacagccagcuuugaugagc
h dp3 AUCUAAGCCUCAUCAAGUGGUUGUGAUAUGGAUACCCAACGCAUAUCACAGCCAGCUUUGAUGAGCUAGGAU
>mir-2b-1
m dm2 uaucacagccagcuuugaggagc
h dm2 CTTCAACTGTCTTCAAAGTGGCAGTGACATGTTGTCAACAATATTCATATCACAGCCAGCTTTGAGGAGCGTTGCGG
m dp3 uaucacagccagcuuugaggagc
h dp3 CUGCGACGCUCUUUAAAGUGGCGGUGACGUGUUGGUAAUAAUAUUCAUAUCACAGCCAGCUUUGAGGAGCGUUGCGG
```

### .scr
The score sheet file is simply a pprint'ed list of dictionaries. All entries
(dictionaries) include the keys 'name', 'totscore' and one or more 'loc_XX'
where XX are strings representing organisms. Additional keys are the keys of the
fdict function dictionary specified in the criteria file.
```
[{'criteria_name_1': 0.654, 'criteria_name_2': -0.5, ... 'loc_dm2': 44, 'loc_dp3': 12, 'name': 'mir-2a-2', 'totscore': 12.974999999999998},
 {'criteria_name_1': 0.44, 'criteria_name_2': -0.2, ... 'loc_dm2': 10, 'loc_dp3': 12, 'name': 'mir-2b-1', 'totscore': 8.472999999999997},
```
This file is easily parsed by direct eval'ing in python:
```
with open(scoresfile, 'r') as input:
    scores = eval(input.read())
```

## MiRscan3
The most recent unforked version of the code, along with documentation, can be
found at:
http://bartellab.wi.mit.edu/softwareDocs/MiRscan3/Introduction.html

- [Lim et al. 2003a] doi:10.1101/gad.1074403
- [Lim et al. 2003b] doi:10.1126/science.1080372
- [Ruby et al. 2007] doi:10.1101/gr.6597907
