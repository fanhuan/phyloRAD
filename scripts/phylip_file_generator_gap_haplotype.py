#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  phylip_file_generator.py
#  This script takes a directory with files seperated by split_library_simrlls.py and merge them into a phylip file for phylogeny reconstruction
#  Copyright 2016 Huan Fan
#  <hfan22@wisc.edu>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  

import sys, gzip, bz2, os, time
import multiprocessing as mp
from optparse import OptionParser

def smartopen(filename,*args,**kwargs):
    if filename.endswith('gz'):
        return gzip.open(filename,*args,**kwargs)
    elif filename.endswith('bz2'):
        return bz2.BZ2File(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)

usage = "usage: %prog [options]"
version = '%prog 20160913.1'
parser = OptionParser(usage = usage, version = version)
parser.add_option("-d", dest = "dataDir", default = 'data',
                  help = "directory containing the data, default = data/")
parser.add_option("-L", dest = "nloci", default = 100, type= int,
                  help = "number of loci to keep, default = 100")
parser.add_option("-n", dest = "weight", default = 1, type= int,
				  help = "weight of a missing loci, default = 1")

(options, args) = parser.parse_args()

outhandle = file(os.path.basename(options.dataDir.rstrip('/'))+'.phylip', 'w')
n = options.weight
###check the data directory:
if not os.path.isdir(options.dataDir):
    print 'Cannot find data directory {}'.format(options.dataDir)
    sys.exit(2)


###Get sample list:
samples = []
for fileName in os.listdir(options.dataDir):
    if os.path.isdir(os.path.join(options.dataDir, fileName)):
        samples.append(fileName)
    else:
        if not fileName.startswith('.'):
            sample = fileName.split(".")[0]
            if sample in samples:
                sample = fileName.split(".")[0]+fileName.split(".")[1]
                if sample in samples:
                    print 'Error, redundant sample or file names. Aborting!'
                    sys.exit(3)
            samples.append(sample)
samples.sort()


#Calculate SBA loci
SBA = {}
for sample in samples:
	filehandle = open(os.path.join(options.dataDir, sample+'.fa'))
	for line in filehandle:
		if line.startswith('>'):
			locus = int(line.split('_')[1].lstrip('locus'))
			if locus < options.nloci:
				if sample in SBA:
					SBA[sample].append(locus)
				else:
					SBA[sample] = [locus]
read_len = len(line.rstrip())
sba_list = list(reduce(set.intersection,map(set,SBA.values())))
sba_not = options.nloci - len(sba_list)
print(sba_list)
#Write phylip file

outhandle.write('%d\t%d\n'%(len(samples),read_len*len(sba_list)+sba_not*n))
for sample in samples:
	if len(sample) < 10:
		outhandle.write(sample+' '*(10-len(sample)))
	else:
		outhandle.write(sample[:9]+' ')
	filehandle = open(os.path.join(options.dataDir, sample+'.fa'))
	line = filehandle.readline()
	i = 0
	while line:
		if line.startswith('>'):
			locus = int(line.split('_')[1].lstrip('locus'))
			if locus < options.nloci:
				if locus in sba_list:
					line = filehandle.readline().rstrip()
					if locus > i:
						gap = locus-i
						outhandle.write('a-'*n*gap)
						i = locus
					outhandle.write(line)
				else:
					print(sample,locus,i)
					if locus > i:
						gap = locus-i
						outhandle.write('a-'*n*gap)
						i = locus
					outhandle.write('t-'*n)
				i += 1
			else:
				break
		line = filehandle.readline()
	if i < options.nloci:
		gap = options.nloci-i
		outhandle.write('a-'*n*gap)
	outhandle.write('\n')
	filehandle.close()
outhandle.close()

		


