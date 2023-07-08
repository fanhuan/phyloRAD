These scripts were used for the simulation part of the paper and not part of the package.
========

### phylip\_file_generator.py
takes a directory with files split by split_library_simrlls.py and merge them into a phylip format alignment for phylogeny reconstruction.

### phylip\_file\_generator_SBAonly.py
identifies SBA loci and write out alignment for them.

### phylip\_file\_generator_gap.py
codes loci that are not shared by all as t---(present) or a---(absent) so those loci are considered by dnadist from phylip.
This is achieved by first identify SBA loci and write out alignment for them. Then t---/a--- for loci that are not SBA.

### phylip\_file\_generator\_gap_V2.py
keeps the sequences at the loci that are not shared by all, while adding t---(present) or a---(absent) so those loci present/absent info are considered by dnadist from phylip.
The difference between this one and v1 is that v1 got rid of the sequences from loci that are SBA completely.

### phylip\_file\_generator_consense.py
takes a directory with files split by split_library_simrlls.py and merge them into a phylip format alignment for phylogeny reconstruction.
  This script is different from phylip_file_generator.py in the sense that it does the consense among the reads from the same loci and same species first, before concatinating them. Any heterozygote site is represented by N.
### phylip\_file\_generator\_consense_N.py
This is a special version of phylip_file_generator_consense.py where it calcuates the proportion of Ns in the alignment.

### split\_libraries\_fastq_FH.py
is used to split Laura's metagenomic data into seperate files for AAF

### split\_libraries\_fastq\_simrlls.py
is used to split fastq files generated from simrlls.

### split\_libraries\_fastq\_simrlls_v2.py
Comparing to the original version:

+ Does not generate a folder containing the original SBA loci (splited for each species)
+ Drop out on loci instead of reads
+ Does not support haplotype

### split\_libraries\_fastq\_simrlls_v3.py
keeps stats of what was there, what’s left after random dropout and what has been selected after reads selection
