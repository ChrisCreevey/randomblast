# randomblast
A tool for identifying gene famililes from very large datasets.


The Randomblast algorithm works by choosing a random sequence from the database and blasting it against all the other sequences in the same database.  
After that, all the sequences with a significant hit are removed from the database and written to a file representing the gene family. 

Another new sequence is randomly picked from the remaining sequences in the database and blasted, and the process is repeated until all sequences are removed from the database. 

This method has previously been shown to work well to define sets of orthologs for supertree reconstruction (Pisani et al 2007), as it efficiently breaks multi-gene families into its paralogy groups.  

Curently, randomblast assumes that youe have provided an fasta formatted file of amino-acid sequences called "dbaa" in the directory from which it is called (already formatted for blast). 

To run the analysis, make sure that you have BLAST (currently randomblast does not support BLAST+, so please use version 2.2.26 - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/legacy/2.2.26/ )
Call randomblast using the command

```
randomblast
```

The files with the sequences belonging to predicted gene famililes will be written to the current directory.

References:


Pisani D, Cotton JA, McInerney JO. Supertrees disentangle the chimerical origin of eukaryotic genomes. Molecular biology and evolution. 2007;24(8):1752-60

