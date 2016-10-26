#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  aaf_distance.py
#  
#  Copyright 2013, 2014,2015 Huan Fan <hfan22@wisc.edu> & Yann Surget-Groba
#  <yann@xtbg.org.cn>
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

import sys, os, math, gzip, time
import multiprocessing as mp
from optparse import OptionParser
from AAF import smartopen, is_exe

def countShared(lines, sn): #count nshare only, for shared kmer table
    shared = [[0] * sn for i in range(sn)]
    for line in lines:
        line = line.split()
        if len(line) == sn+1:
            line = line[1:]
        line = [int(i) for i in line]
        for i in range(sn):
            for j in range(i + 1, sn):
                if line[i]*line[j] != 0:
                    shared[i][j] += 1
    return shared


def aaf_dist(datfile,countfile,nThreads,samples,kl,long=False,)
    #check executables
	if os.system('which fitch_kmerX > /dev/null'):
		if long:
			fitch = './fitch_kmerX_long'
		else:
			fitch = './fitch_kmerX'
		if not is_exe(fitch):
			print(fitch+' not found. Make sure it is in your PATH or the')
			print('current directory, and that it is executable')
			sys.exit()
	else:
		if long:
			fitch = 'fitch_kmerX_long'
		else:
			fitch = 'fitch_kmerX'
	#process the .dat.gz file
	try:
		iptf = smartopen(datfile,'rt')
	except IOError:
		print('Cannot open file', datafile)
		sys.exit()

	if not os.path.isfile(countfile):
		print('Cannot find file', countfile)
		sys.exit()

	try:
		total = open(countfile,'rt')
	except IOError:
		print('Cannot open file', countfile)
		sys.exit()

	try:
		infile = open('infile','wt')
	except IOError:
		print('Cannot open infile for writing')
		sys.exit()

	###Initialize shared kmers matrix
	sn = len(samples)    #species number
	nshare = [[0] * sn for i in range(sn)]

	###Compute the number of lines to process per thread
	line = iptf.readline()
	line_size = sys.getsizeof(line)
	chunkLength = int(1024 ** 3 / nThreads / line_size)
	print('chunkLength = {}'.format(chunkLength))

	###Compute shared kmer matrix
	nJobs = 0
	pool = mp.Pool(nThreads)
	results = []
	print('{} start running jobs'.format(time.strftime('%c')))
	print('{} running {} jobs'.format(time.strftime('%c'), nThreads))
	while True:
		if nJobs == nThreads:
			pool.close()
			pool.join()
			for job in results:
				shared = job.get()
				for i in range(sn):
					for j in range(i + 1, sn):
						nshare[i][j] += shared[i][j]
			pool = mp.Pool(nThreads)
			nJobs = 0
			results = []
			print('{} running {} jobs'.format(time.strftime('%c'), nThreads))

		lines = []
		for nLines in range(chunkLength):
			if not line: #if empty
				break
			lines.append(line)
			line = iptf.readline()
		if not lines: #if empty
			break 
		job = pool.apply_async(countShared, args=[lines, sn])
		
		results.append(job)
		nJobs += 1

	if nJobs:
		print('{} running last {} jobs'.format(time.strftime('%c'), len(results)))
		pool.close()
		pool.join()
		for job in results:
			shared = job.get()
			for i in range(sn):
				for j in range(i + 1, sn):
					nshare[i][j] += shared[i][j]

	iptf.close()

	###Compute distance matrix
	ntotal = [0.0] * sn

	for i in range(sn):
		ntotal[i] = float(total.readline().split()[1])
	dist = [[0] * sn for i in range(sn)]    

	for i in range(sn):
		for j in range(i + 1, sn):
			mintotal = min(ntotal[i], ntotal[j])
			if nshare[i][j] == 0:
				dist[j][i] = dist[i][j] = 1
			else:
				distance = (-1 / float(kl) * math.log(nshare[i][j] / mintotal)
				dist[j][i] = dist[i][j] = distance
				nshare[j][i] = nshare[i][j]

	total.close()

	###Write infile
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
	fh1 = open('aaf.tre','wt')

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
	command = 'mv infile aaf.dist'
	os.system(command)

	os.system('rm -f outfile outtree')

	print('{} end'.format(time.strftime("%c")))
