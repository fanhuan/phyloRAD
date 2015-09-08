#!/usr/bin/env python

import sys, re, subprocess, bisect, os, collections, random
from Bio import SeqIO
from optparse import OptionParser
from scipy.stats import poisson

#Get command-line arguments
Usage = "%prog [options] -i <input filename>"
version = '%prog 20150903.1'
parser = OptionParser()
parser.add_option('-i', dest='infile', help='prefix of infiles', default=False)
parser.add_option('-d', dest='outdir', help='directory for output files',
                  default=False)
parser.add_option('-e', dest='enzyme', help='restriction enzyme, SbfI(default), EcoRI, PstI, ApeKI or MseI',
                  default='SbfI')
parser.add_option('-r', dest='rate', type = float, help='rate of drop one tag',
                  default=0.1)
parser.add_option('-a', dest='all', action='store_true',
                  help='out put shared-by-all sites only', default=False)
parser.add_option('-n', dest='noOut', action='store_true', help='only count sites',
                  default=False)
parser.add_option('-v', dest='verbose', action='store_true', help='verbose',
                  default=False)
(options, args) = parser.parse_args()

infile   = options.infile
outdir  = options.outdir
noOut    = options.noOut
verbose = options.verbose

#Define restriction enzyme recognition site

SiteDic = {'SbfI':['CCTGCAGG'],'EcoRI':['GAATTC'],'PstI':['CTGCAG'],
           'ApeKI':['GCAGC','GCTGC'],'MseI':['TTAA']}
cutSite = SiteDic[options.enzyme]


#parse genome fasta file to get all rad-tags
if verbose:
    print 'Reading input sequences'

if outdir == False:
    print('Quitting the program. No output directory specified.')
    sys.exit(2)


#output directory

if outdir in os.listdir('./'):
    if verbose:
        answer = raw_input('{} exists, overwrite? Y/N\n'.format(outdir))
        if answer == 'Y' or answer =='y':
            print('{} is going to be overwritten'.format(outdir))
            os.system('rm -r -f {}/*'.format(outdir))
            os.system('rm -f {}/*'.format(outdir))
        elif answer == 'N' or answer =='n':
            print('Quitting the program. Please rerun with new output directory.')
            sys.exit(2)
        else:
            print('Wrong keyboard input, exit')
            sys.exit(2)
    else:
        os.system('rm -r -f {}/*'.format(outdir))
        os.system('rm -f {}/*'.format(outdir))
else:
    os.system('mkdir {}'.format(outdir))

alignment = open(infile+'.phy')

sample = [] #sample list
dic = {} #dash dictionary dictionary
dashtotal = {}

if verbose:
    print 'construct the dash dictionary: where the dashes are'

# read the number of species on the second line
sn = int(alignment.readline().split()[0])
lines = alignment.readlines()
block = (len(lines)-1)/(sn+1)
alignment.close()
alignment = open(infile+'.phy')
# read the first sn rows with sample names
alignment.readline() #skip the first row with the number of samples and characters
for i in range(sn):
    row = alignment.readline().split()
    sample.append(row[0])
    dash = [match.start() for match in re.finditer('-',row[1])]
    if dash != []: # if there are dashed in the row
        dic[sample[i]] = {dash[0]+1:len(dash)}
        dashtotal[sample[i]] = len(dash)
    else:
        dic[sample[i]] = {}
        dashtotal[sample[i]] = 0

# read the rest of the alignment, each block contains sn+1(blank line) lines
for j in range(1,block):
    alignment.readline() #skip the empty line
    for i in range(sn):
        row = alignment.readline().lstrip()
        if row != '': #check if we really skipped the empty line
            dash = [match.start() for match in re.finditer('-',row)]
            if dash != []:
                dic[sample[i]][dash[0]+ 1 + 60*j- dashtotal[sample[i]]]= dashtotal[sample[i]]+len(dash)
                dashtotal[sample[i]] = dashtotal[sample[i]]+len(dash)
        else:
            if verbose:
                print 'Did not skip the empty line properly at block'
                print j,i
alignment.close()

if verbose:
    print 'Species list:\n'
    print sample
    print 'START PROCESSING THE FASTA FILE'
handle = open(infile+'.fas')
dashDic_r = {} # real
dashDic_s = {} # after drop some random site
for record in SeqIO.parse(handle, 'fasta'):
    tags = []
    s = str(record.seq)
    s = s.upper()
    sitelist = sorted(dic[record.id].keys())
    dashlist = [] #restriction sites in alignment
    selectedlist = [] #seq ID of output reads (not dropped)
    #cut current sequence at all recognition sites
    fragments = re.split(cutSite[0], s)
    matches = re.finditer(cutSite[0],s)
    matchSites = [] #restriction sites in fasta
    if len(cutSite) > 1: #ApeKI
        for cutsite in cutSite:
            fragments += re.split(cutSite[1], s)
            matches += re.finditer(cutSite[1],s)
    for match in matches:
        matchSites.append(match.start())
    matchSites = sorted(matchSites)
    for matchsite in matchSites:
        if len(dic[record.id]) == 0: #no dashes in the alignment
            dashlist.append(matchsite)
        else:
            dashlist.append(matchsite + dic[record.id][sitelist[bisect.bisect(sitelist,matchsite)-1]])
    nCut = len(fragments) - 1
    #if there are recognition sites
    taglist = []
    if len(fragments) > 1:
        first = fragments[0]
        last  = fragments[-1]
        if len(first) > 100:
            taglist.append(str(dashlist[0])+'_L')
            tags.append(first[-100:])
        
        #we process internal fragments (restriction site on each side)
        fragments = fragments[1:-1]
        k = 0
        for f in fragments:
            if len(f) > 100:
                tags.append(f[:100])
                taglist.append(str(dashlist[k])+'_R')
                tags.append(f[-100:])
                k += 1
                taglist.append(str(dashlist[k])+'_L')
            #if the fragment is too short, we still need to move to the next
            #restriction site
            else:
                k += 1
        
        
        #we take the first 100bp of the last fragment if it is long enough
        if len(last) > 100:
            taglist.append(str(dashlist[-1])+'_R')
            tags.append(last[:100])
    if verbose:
        print record.id
        print 'Number of cut sites:', nCut
        print 'Number of rad-tags:', len(tags)

    #write reads to file
    if not noOut:
        if verbose:
            print 'Write tags to file'
        output = open(outdir+'/'+record.id+'_RAD.fa', 'w')
        m = 0
        for read in tags:
            if options.rate < random.random():
                output.write('>%s%s%s\n%s\n' % (record.id, '_', taglist[m], read))
                selectedlist.append(taglist[m])
            m += 1
        output.close()
    dashDic_r[record.id] = dashlist
    dashDic_s[record.id] = selectedlist

handle.close()
# Write the 0/1 of restriction sites to a file
rest = open(infile + '_rest.txt','w')
rest.write(' ')
select = open(infile + '_select.txt','w')
select.write(' ')

fullset = sorted(set().union([item for sublist in dashDic_r.values() for item in sublist]))

selectset = sorted(set().union([item for sublist in dashDic_s.values() for item in sublist]))

for taxa in sample:
    rest.write('\t%s'%(taxa))
    select.write('\t%s'%(taxa))
rest.write('\t%s\n'%('Total'))
select.write('\t%s\n'%('Total'))
sum_r = []
sum_s = []


for item in fullset:
    rest.write('%s'%(item))
    total = 0
    for taxa in sample:
        if item in dashDic_r[taxa]:
            rest.write('\t%d'%(1))
            total += 1
        else:
            rest.write('\t%d'%(0))
    sum_r.append(total)
    rest.write('\t%d\n'%total)

summary = {}
for item in selectset:
    select.write('%s'%(item))
    total = 0
    for taxa in sample:
        if item in dashDic_s[taxa]:
            select.write('\t%d'%(1))
            total += 1
        else:
            select.write('\t%d'%(0))
    sum_s.append(total)
    select.write('\t%d\n'%total)
    summary[item] = total


hist1 = collections.Counter(sum_r)
rest.close()
hist2 = collections.Counter(sum_s)
select.close()

if verbose:
    print 'In the alignment:'
    print hist
    print float(hist1[sn])/len(fullset)*100,' percent of restriction sites are shared by all.'


    print 'In the output reads files:'

    print hist
    print float(hist2[sn])/len(selectset)*100,' percent of restriction sites are shared by all.'

else:
    print float(hist1[sn])/len(fullset)*100
    print float(hist2[sn])/len(selectset)*100

if options.all:
    os.system('mkdir {}_all'.format(outdir))
    for id in sample:
        handle_in = open(outdir+'/'+id+'_RAD.fa')
        handle_out = open(outdir+'_all/'+id+'_RAD.fa','w')
        for record in SeqIO.parse(handle_in, 'fasta'):
            if summary[record.id[len(id)+1:]] == len(sample):
                handle_out.write('>%s\n%s\n' % (record.id, record.seq))