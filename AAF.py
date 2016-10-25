#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  aaf_pairwise.py
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

import sys, gzip, bz2, os, time
import multiprocessing as mp

def smartopen(filename,*args,**kwargs):
    if filename.endswith('gz'):
        return gzip.open(filename,*args,**kwargs)
    elif filename.endswith('bz2'):
        return bz2.BZ2File(filename,*args,**kwargs)
    else:
        return open(filename,*args,**kwargs)

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def aaf_kmercount(dataDir,samples,k,n,nThreads,memPerThread):
    #check excutables
    if k > 25:
        if os.system('which kmer_countx > /dev/null'):
            kmerCount = './kmer_countx'
            if not is_exe(kmerCount):
                print('kmer_countx not found. Make sure it is in your PATH or the')
                print('current directory, and that it is executable')
                sys.exit(1)
        else:
            kmerCount = 'kmer_countx'

    else:
        if os.system('which kmer_count > /dev/null'):
            kmerCount = './kmer_count'
            if not is_exe(kmerCount):
                print('kmer_count not found. Make sure it is in your PATH or the')
                print('current directory, and that it is executable')
                sys.exit(1)
        else:
            kmerCount = 'kmer_count'

    ###Prepare kmer_count jobs
    jobList = []
    for sample in samples:
        outFile = '{}.pkdat.gz'.format(sample)
        command = '{} -l {} -n {} -G {} -o {} -f '.format(kmerCount, k, n, memPerThread, outFile)
        command1 = ''
        for inputFile in os.listdir(os.path.join(dataDir, sample)):
            inputFile = os.path.join(dataDir, sample, inputFile)
            handle = smartopen(inputFile)
            firstChar = handle.read(1)
            if firstChar == '@':
                seqFormat = 'FQ'
            elif firstChar == '>':
                seqFormat = 'FA'
            else:
                print('Error, file {} is not FA or FQ format. Aborting!'.format(inputFile))
                sys.exit(3)
            command1 += " -i '{}'".format(inputFile)
            command += '{}{}'.format(seqFormat,command1)
            jobList.append(command)
    jobList = jobList[::-1]

    ###Run jobs
    pool = mp.Pool(nThreads)
    jobs = []
    nJobs = 0
    batch = 0
    count = 0
    nBatches = int(len(jobList) / nThreads)
    if len(jobList) % nThreads:
        nBatches += 1

    while 1:
        if nJobs == nThreads:
            batch += 1
            print(time.strftime('%c'))
            print("running batch {}/{}".format(batch, nBatches))
            for job in jobs:
                print(command)
                pool.apply_async(os.system(command))
            pool.close()
            pool.join()
            pool = mp.Pool(nThreads)
            nJobs = 0
            jobs = []
        if jobList:
            command = jobList.pop()
            jobs.append(command)
            #job = pool.apply_async(runJob, args=[command, options.sim])
            nJobs += 1
        else:
            break
        count += 1

    if nJobs:
        print(time.strftime('%c'))
        print("running last batch")
        for job in jobs:
            pool.apply_async(os.system(command))
        pool.close()
        pool.join()


