#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  phyloRAD_pairwise.py
#
#  Copyright 2016 Huan Fan <hfan22@wisc.edu>
#
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

import sys, gzip, bz2, os, time, subprocess, math
import multiprocessing as mp
from optparse import OptionParser
from AAF import *


usage = "usage: %prog [options]"
version = '%prog 20170913.1'
parser = OptionParser(usage = usage, version = version)
parser.add_option("-k", dest = "kLen", type = int, default = 25,
                  help = "k for reconstruction, default = 25")
parser.add_option("--ks", dest = "ksLen", type = int, default = 25,
                  help = "k for reads selection, default = 25")
parser.add_option("-n", dest = "filter", type = int, default = 1,
                  help = "k-mer filtering threshold, default = 1")
parser.add_option("-d", dest = "dataDir", default = 'data',
                  help = "directory containing the data, default = data/")
parser.add_option("-G", dest = "memSize", type = int, default = 4,
                  help = "total memory limit (in GB), default = 4")
parser.add_option("-t", dest = "nThreads", type = int, default = 1,
                  help = "number of threads to use, default = 1")
parser.add_option("-l", dest = "long", action = 'store_true',
                  help = "use fitch_kmerX_long instead of fitch_kmerX")

(options, args) = parser.parse_args()

n = options.filter
memSize = options.memSize
nThreads = options.nThreads
kl = options.kLen
ks = options.ksLen
memPerThread = int(options.memSize / float(nThreads))
dataDir = options.dataDir

if not memPerThread:
	print('Not enough memory, decrease nThreads or increase memSize')
	sys.exit()

###check the data directory:
if not os.path.isdir(dataDir):
    print('Cannot find data directory {}'.format(dataDir))
    sys.exit(2)


###check for the executable files:
#kmer_merge
if os.system('which kmer_merge > /dev/null'):
    filt = './kmer_merge'
    if not is_exe(filt):
        print('kmer_merge not found. Make sure it is in your PATH or the')
        print('current directory, and that it is executable')
        sys.exit(1)
else:
    filt = 'kmer_merge'

#ReadsSelector
if os.system('which ReadsSelector > /dev/null'):
    ReadsSelector = './ReadsSelector'
    if not is_exe(filt):
        print('ReadsSelector not found. Make sure it is in your PATH or the')
        print('current directory, and that it is executable')
        sys.exit(1)
else:
    ReadsSelector = 'ReadsSelector'

#fitch
if os.system('which fitch_kmerX > /dev/null'):
    if options.long:
        fitch = './fitch_kmerX_long'
    else:
        fitch = './fitch_kmerX'
    if not is_exe(fitch):
        print(fitch+' not found. Make sure it is in your PATH or the')
        print('current directory, and that it is executable')
        sys.exit()
else:
    if options.long:
        fitch = 'fitch_kmerX_long'
    else:
        fitch = 'fitch_kmerX'


#set up directory for selected reads
selection_dir = '{}_ks{}_pairwise'.format(os.path.basename(dataDir.rstrip('/')),ks)

if os.path.exists('./'+selection_dir):
	command = 'rm -r {}'.format(selection_dir)
	os.system(command)
command = 'mkdir {}'.format(selection_dir)
os.system(command)


###Run aaf_kmercount to get pkdat for each species

samples = aaf_kmercount(dataDir,ks,n,nThreads,memPerThread)

###Build distance matrix
dist = [[0] * sn for i in range(sn)]
for i in range(sn):
    for j in range(i+1,sn):
        #Reads selection
        command = []
        command.append('mkdir {}/{}_{}'.format(selection_dir,samples[i],samples[j]))
        #command.append = 'mkdir {}/{}_{}/{}'.format(selection_dir,samples[i],samples[j],samples[i])
        #command.append = 'mkdir {}/{}_{}/{}'.format(selection_dir,samples[i],samples[j],samples[j])
        command.append('{} -k s -c -d 0 -A A -B A {}.pkdat.gz {}.pkdat.gz | cut -f 1 > test.kmer'.format(filt,samples[i],samples[j]))
        command.append('{} -k test.kmer -fa 1 -o {}/{}_{}/{} -s {}/{}/*' \
                  .format(ReadsSelector,selection_dir,samples[i],samples[j],
                          samples[i],dataDir,samples[i]))
        command.append('{} -k test.kmer -fa 1 -o {}/{}_{}/{} -s {}/{}/*' \
            .format(ReadsSelector,selection_dir,samples[i],samples[j],
                    samples[j],dataDir,samples[j]))
        for comm in command:
            print(comm)
            os.system(comm)
        #kmer_count
        ntotal = []
        command = '{} -l {} -n {} -G {} -o {}_temp.pkdat -f FA -i {}/{}_{}/{}.fa'.format(kmerCount,kl,n,memSize/2,samples[i],selection_dir,samples[i],samples[j],samples[i])
        output = subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT)
        ntotal.append(float(output.decode('ascii').split()[1]))
        command = '{} -l {} -n {} -G {} -o {}_temp.pkdat -f FA -i {}/{}_{}/{}.fa'.format(kmerCount,kl,n,memSize/2,samples[j],selection_dir,samples[i],samples[j],samples[j])
        output = subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT)
        ntotal.append(float(output.decode('ascii').split()[1]))
        #kmer_merge
        command = "{} -k s -c -d '0' -a 'T,M,F' {}_temp.pkdat {}_temp.pkdat | wc -l".format(filt,samples[i],samples[j])
        output = subprocess.check_output(command,shell=True,stderr=subprocess.STDOUT)
        nshared = int(output.decode('ascii').split()[0])
        if nshared == 0:
            distance = 1
        else:
            distance = (-1.0 / kl) * math.log(nshared / min(ntotal))
        dist[j][i] = dist[i][j] = distance

os.system('rm *.pkdat*')

###constructing the tree
###Write infile
try:
	infile = open('infile','w')
except IOError:
	print('Cannot open infile for writing')
	sys.exit()

infile.write('{} {}'.format(sn, sn))
namedic = {}
for i in range(sn):
    lsl = len(samples[i])
    if lsl >= 10:
        ssl = samples[i][:10]
        appendix = 1
        while ssl in namedic:
            if appendix < 10:
                ssl = samples[i][:9]+str(appendix)
            elif appendix > 9:
                ssl = samples[i][:8]+str(appendix)
            appendix += 1
    else:
        ssl = samples[i] + ' ' * (10 - lsl)
    namedic[ssl] = samples[i]
    infile.write('\n{}'.format(ssl))
    for j in range(sn):
        infile.write('\t{}'.format(dist[i][j]))

infile.close()

###Run fitch_kmer
print('{} building tree'.format(time.strftime("%c")))
if os.path.exists("./outfile"):
    os.system("rm -f outfile outtree")
command = 'printf "K\n{}\nY" | {} > /dev/null'.format(int(kl),fitch)
os.system(command)
fh = open('outtree','rt')
fh1 = open(selection_dir +'.tre','wt')

for line in fh:
    for key in namedic:
        key_new = key.rstrip()+":"
        if key_new in line:
            newline = line.replace(key_new,namedic[key].rstrip()+":",1)
            line = newline
    fh1.write(line) #This can be either line or new line because when it exits
                    #the for loop, line==newline
fh.close()
fh1.close()
command = 'mv infile {}.dist'.format(selection_dir)
os.system(command)

os.system('rm -f outfile outtree')

print('{} end'.format(time.strftime("%c")))
