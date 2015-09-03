#!/usr/bin/env python

import sys, re, subprocess, bisect, os, collections
from Bio import SeqIO
from optparse import OptionParser
from scipy.stats import poisson

#Get command-line arguments
Usage = "%prog [options] -i <input filename>"
version = '%prog 20150803.1'
parser = OptionParser()
parser.add_option('-d', dest='dataDir', help='directory with the selected reads', default=False)
(options, args) = parser.parse_args()

###Get sample list:
samples = []
sites = {}
for fileName in os.listdir(options.dataDir):
    sample = fileName.split('.')[0]
    samples.append(sample)
    handle = open(os.path.join(options.dataDir, fileName))
    sites[sample] = []
    for record in SeqIO.parse(handle, 'fasta'):
        sites[sample].append('_'.join(record.id.split('_')[-2:]))
    handle.close()
# Write the 0/1 of restriction sites to a file
rest = open(options.dataDir + '_rest.txt','w')
rest.write(' ')

fullset = sorted(set().union([item for sublist in sites.values() for item in sublist]))

samples.sort()

for taxa in samples:
    rest.write('\t%s'%(taxa))
rest.write('\t%s\n'%('Total'))

sum = []

for item in fullset:
    rest.write('%s'%(item))
    total = 0
    for taxa in samples:
        if item in sites[taxa]:
            rest.write('\t%d'%(1))
            total += 1
        else:
            rest.write('\t%d'%(0))
    sum.append(total)
    rest.write('\t%d\n'%total)

hist = collections.Counter(sum)
print hist
print len(fullset)
print float(hist[len(samples)])/len(fullset)*100,' percent of restriction sites are shared by all.'
rest.close()


