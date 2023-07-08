**phyloRAD**
===========

A package reconstructs phylogeny from RADseq data using alignment and assembly-free (AAF) method. 

__Citation__: 
Reconstructing phylogeny from reduced‐representation genome sequencing data without assembly or alignment, Huan Fan, Anthony R Ives, Yann Surget‐Groba, June 2018, Molecular Ecology Resources
DOI: 10.1111/1755-0998.12921, [Full text](https://rdcu.be/6cok)


This package includes executables from Phylokmer (Author: Jue Ruan), PHYLIP (Author: Joe Felsenstein) and ReadsSelector (Author: Chengxi Ye).

Requirements
------------
+ UNIX system (Linux or macOS)
+ g++/gcc compilers
+ Python 2.7+ or higher (including 3.x)

__Note that only v1.0 supports both python2.7 and python 3.x. v1.1+ only supports python3.5 and 3.6.__


Installation:
-------------
Warning: the new version (where gzip files are taken) is only tested on linux.  
1. Download or clone this repo.  
2. Install a c++ library called gzstream so that ReadsSelector could process gzipped files.  
	i. `cd` into the subdirectory `source`   
	ii. `tar xzvf gzstream.tgz`  
	iii. `cd gzstream`
	iv. Run `make`  
	v. Run `export CPLUS_INCLUDE_PATH=PATH_TO_phyRAD/source/gzstream`  
	vi. Run `export LIBRARY_PATH=PATH_TO_phyRAD/source/gzstream`
3. `cd ../` out into `source`, run `make all` and `make clean`  
4. `cp` the executables you just compiled to your PATH or working directory, including consense, treedist, fitch_kmerX, fitch_kmerX_long, kmer_count(x), kmer_merge and ReadsSelector. 
5. Make sure python knows where the AAF module is. You can achive this by 
	+ Put the AAF.py in your working directory
	+ Add where AAF.py is to your PATH `export PATH=/PATH/TO/AAF.py:$PATH`
	+ Put AAF.py in your site-packages folder. If you don't know where your site-packages folder is, try: `>>> import site; site.getsitepackages()`

Usage and Examples: 
---------------
__Input__: Fasta or fastq files. If there are multiple files for one sample, put them in one folder.  
__Output__: Newick format phylogenetic tree.

Using __shared-by-all__ reads selection::
 

    $ phyloRAD_sba.py -h
    
    Usage: phyloRAD_sba.py [options]

	Options:
	  --version    show program's version number and exit
	  -h, --help   show this help message and exit
	  -k KLEN      k for reconstruction, default = 25
	  --ks=KSLEN   k for reads selection, default = 25
	  -n FILTER    k-mer filtering threshold, default = 1
	  -d DATADIR   directory containing the data, default = data/
	  -G MEMSIZE   total memory limit (in GB), default = 4
	  -t NTHREADS  number of threads to use, default = 1
	  -l           use fitch_kmerX_long instead of fitch_kmerX

__Example__:
	
	 $ phyloRAD_sba.py -d test_data -G 4 -t 2 -k 15 --ks 15 -n 2
	 
Using __pairwise__ reads selection::

    $ phyloRAD_pairwise.py -h
	Usage: phyloRAD_pairwise.py [options]
	
	Options:
	  --version    show program's version number and exit
	  -h, --help   show this help message and exit
	  -k KLEN      k for reconstruction, default = 25
	  --ks=KSLEN   k for reads selection, default = 25
	  -n FILTER    k-mer filtering threshold, default = 1
	  -d DATADIR   directory containing the data, default = data/
	  -G MEMSIZE   total memory limit (in GB), default = 4
	  -t NTHREADS  number of threads to use, default = 1
	  -l           use fitch_kmerX_long instead of fitch_kmerX
 

__Examples__:

	$ phyloRAD_pairwise.py -d test_data -G 4 -t 2 -k 15 --ks 15 -n 2
