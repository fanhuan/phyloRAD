**phyloRAD**
===========

A package reconstructs phylogeny from RADseq data using alignment and assembly-free (AAF) method.

This package includes executables from Phylokmer (Author: Jue Ruan), PHYLIP (Author: Joe Felsenstein) and ReadsSelector (Author: Chengxi Ye).

Requirements
------------
+ UNIX system (Linux or macOS)
+ g++/gcc compilers
+ Python 2.7+ or higher (including 3.x)


Installation:
-------------
1. Download or clone this repo.  
2. 'cd' into the top-level `source/` directory 
3. Run 'make all'  and 'make clean'
4. 'cp' the executables you just compiled to your PATH or working directory, including consense, treedist, fitch_kmerX, fitch_kmerX_long, kmer_count(x), kmer_merge and ReadsSelector. 
5. Make sure python knows where the AAF module is. You can achive this by 
	+ Put the AAF.py in your working directory
	+ Add where AAF.py is to your PATH `export PATH=/PATH/TO/AAF.py:$PATH`
	+ Put AAF.py in your site-packages folder. If you don't know where your site-packages folder is, try: `>>> import site; site.getsitepackages()`

Usage and Examples: 
---------------

Using shared-by-all reads selection::
 

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


Using pairwise reads selection::

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
 

Examples:

The command used in generating Figure 5c in the phyloRAD paper:


    $ phyloRAD_sba.py -d Quercus_CAVENDER-BARES -G 20 -t 5 -k 21 --ks 21 -n 2
    
The command used in generating Figure 5d in the phyloRAD paper:  

	$ phyloRAD_pairwise.py -d Quercus_CAVENDER-BARES -G 20 -t 5 -k 21 --ks 21 -n 2
