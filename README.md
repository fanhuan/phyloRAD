**AAF-RAD**
===========

A package reconstructs phylogeny from RADseq data using alignment and assembly-free (AAF) method.

This package includes executables from Phylokmer (Author: Jue Ruan) and PHYLIP (Author: Joe Felsenstein) and ReadsSelector (Author: Chengxi Ye).

Requirements
------------
+ UNIX system (Linux or macOS)
+ g++/gcc compilers
+ Python 2.7+


Installation:
-------------
1. Download or clone this repo.  
2. 'cd' into the top-level `source/` directory 
3. Run 'make all'  and 'make clean'
4. 'cp' the executables you just compiled to your PATH or working directory, including consense, treedist, fitch_kmerX, fitch_kmerX_long, kmer_count(x), kmer_merge and ReadsSelector. 

Example usage: 
---------------

See all parameter options:

.. code:: bash  

    $ AAF-RAD_sba.py -h
    
    optional arguments:
    -h, --help      show this help message and exit
    --version       show program's version number and exit
    -o outname      [str] output file name prefix (default 'out')


Modified parameters::

    $ AAF-RAD_sba.py -o test2 -ks 15 

.. code:: bash  

    $ echo "hellow world" > treefile   





0
