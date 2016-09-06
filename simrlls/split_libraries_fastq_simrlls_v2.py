#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  split_libraries_fastq_FH.py
#
#  This script is written for spliting simrlls simulation data. This also includes 1) getting rid of barcode 2) introduce random dropout rate
#  
#  Copyright 2016 Huan Fan <hfan22@wisc.edu>
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

import sys, gzip, bz2, os, time, math, re, argparse,random
import multiprocessing as mp
from optparse import OptionParser

def smartopen(filename,*args,**kwargs):
	'''opens with open unless file ends in .gz, then use gzip.open
		in theory should transparently allow reading of files regardless of
		compression'''
	if filename.endswith('gz'):
		return gzip.open(filename,*args,**kwargs)
	elif filename.endswith('bz2'):
		return bz2.BZ2File(filename,*args,**kwargs)
	else:
		return open(filename,*args,**kwargs)

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)



usage = "usage: %prog [args]"
version = '%prog 20160726.1'
parser = argparse.ArgumentParser()
parser.add_argument("-i", dest = "input", action='append')
parser.add_argument("-r", dest = "rate", type = float, default = 0,
                  help = "dropout rate")
parser.add_argument("-d", dest = "dir", default = "test",
                  help = "output directory for fasta files after seperation, default = test")
parser.add_argument("-L", dest = "nloci", type =  int, default = 100,
						help = "number of loci to sample, default = 100")
args = parser.parse_args()

rate = args.rate
outputDir = args.dir

###Set up output folder:
if os.path.exists(outputDir+'_r'+str(rate)):
    print('The output directory {}_r{} already exist.'.format(outputDir,str(rate)))
    print('It is going to be over written.')
    os.system('rm -r {}_r{}'.format(outputDir,str(rate)))

os.system('mkdir {}_r{}'.format(outputDir,str(rate)))

### Start processing input file
# Make the before selection directory
from Bio import SeqIO
samples = {} # {sample name: sample output file handle}
sba = {} #{sample name: list of locus}

if len(args.input) == 1:
	input_handle = smartopen(args.input[0])
	for seq_record in SeqIO.parse(input_handle,"fastq"):
		locus_id = int(seq_record.id.split('_')[1].lstrip('locus'))
		if locus_id >= args.nloci:
			break
		else:
			sample = seq_record.id.split('_')[2]
			if sample in sba:
				sba[sample].append(locus_id)
			else:
				sba[sample]=[locus_id]
	input_handle.close()

if len(args.input) == 2:
	f1 = smartopen(args.input[0],'rt')
	f2 = smartopen(args.input[1],'rt')
	recs1 = SeqIO.parse(f1, 'fastq')
	recs2 = SeqIO.parse(f2, 'fastq')
	for rec1, rec2 in zip(recs1, recs2):
		locus_id = int(rec1.id.split('_')[1].lstrip('locus'))
		if locus_id >= args.nloci:
			break
		else:
			sample = rec1.id.split('_')[2]
			if sample in sba:
				sba[sample].append(locus_id)
			else:
				sba[sample]=[locus_id]
	f1.close()
	f2.close()

# Add random dropout
sba_selected = {}
for sample in sba:
	loci_list = set(sba[sample])
	loci_list = random.sample(loci_list,int(len(loci_list)*(1-rate)))
	sba_selected[sample] = loci_list


# Write fasta files
out_sample = {}
if len(args.input) == 1:
	input_handle = smartopen(args.input[0],'rt')
	for seq_record in SeqIO.parse(input_handle,"fastq"):
		sample = seq_record.id.split('_')[2]
		locus_id = int(seq_record.id.split('_')[1].lstrip('locus'))
		if locus_id >= args.nloci:
			break
		else:
			if locus_id in sba_selected[sample]:
				if sample in out_sample:
						out_sample[sample].write('>'+seq_record.id+'\n')
						out_sample[sample].write(str(seq_record.seq[6:])+'\n')
				else:
					out_sample[sample] = open(outputDir+'_r'+str(rate)+'/'+sample+'.fa','w')
					out_sample[sample].write('>'+seq_record.id+'\n')
					out_sample[sample].write(str(seq_record.seq[6:])+'\n')

if len(args.input) == 2:
	f1 = smartopen(args.input[0],'rt')
	f2 = smartopen(args.input[1],'rt')
	recs1 = SeqIO.parse(f1, 'fastq')
	recs2 = SeqIO.parse(f2, 'fastq')
	for rec1, rec2 in zip(recs1, recs2):
		sample = rec1.id.split('_')[2]
		locus_id = int(rec1.id.split('_')[1].lstrip('locus'))
		if locus_id in sba_selected[sample]:
			if sample+'_1' in out_sample:
				out_sample[sample+'_1'].write('>'+rec1.id+'\n')
				out_sample[sample+'_1'].write(str(rec1.seq[6:])+'\n')
				out_sample[sample+'_2'].write('>'+rec2.id+'\n')
				out_sample[sample+'_2'].write(str(rec2.seq[6:])+'\n')
			else:
				os.system('mkdir {}_r{}/{}'.format(outputDir,str(rate),sample))
				out_sample[sample+'_1'] = open(outputDir+'_r'+str(rate)+'/'+sample+'/'+sample+'_R1.fa','w')
				out_sample[sample+'_1'].write('>'+rec1.id+'\n')
				out_sample[sample+'_1'].write(str(rec1.seq[6:])+'\n')
				out_sample[sample+'_2'] = open(outputDir+'_r'+str(rate)+'/'+sample+'/'+sample+'_R2.fa','w')
				out_sample[sample+'_2'].write('>'+rec2.id+'\n')
				out_sample[sample+'_2'].write(str(rec2.seq[6:])+'\n')
	f1.close()
	f2.close()

for key in out_sample:
    out_sample[key].close()