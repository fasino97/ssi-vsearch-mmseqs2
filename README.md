# ssi-vsearch-mmseqs2
This repository contains all of the files necessary to run accuracy, homogeneity, and completeness using SSI-VSEARCH and SSI-MMseqs2. 

## Databases
The two databases used in the accompanying paper are the Ribosmal Database Project (RDP) and the SILVA rRNA database project:

The RDP database: https://www.glbrc.org/data-and-tools/glbrc-data-sets/ribosomal-database-project.

The SILVA database: https://www.arb-silva.de/.

## Pre-Processing
Before you can run any incremental clustering algoirthms you must first use the SplitData.py file to sort the data into testing and training data. VSEARCH expects the input file to be a FASTA file, but MMseqs2 requires formatteed databases that can be converted using MMSeqs2 using createdb.py.

## SSI-VSEARCH
SSI-VSEARCH is run in Matlab. The main code is run in VSearchtest.m, where you specify the similarity threshold (75-97), taxonomic depth (2-6), associated map, and the vsearch and usearch executable files. Taxonomic depth does not need to be set for running SSI-VSEARCH, as this is only used in the accuracy calculations available when running the IncrementalVSearch.m file. For this improved version of SSI-VSEARCH we use a more specific post-processing code. Talk about creating map and associated utility files.

All of the information regarding VSEARCH can be found here: https://github.com/torognes/vsearch.

USEARCH is used as well in our SSI-VSEARCH implementation. The download and further information can be found here: https://github.com/rcedgar/usearch12.

## SSI-MMSeqs2
SSI-MMSeqs2 is run in Python. This code was run mainly on Rowan University's HPC so the code attached (slurm-script) is the code to call MMSeqs2 and capture the time it takes to run. A database must be used as the input by using MMSeqs2's createdb function, and batches are added using MMSeqs2's clusterupdate function.

The download and inforamtion about MMSeqs2 can be found on this Github repository: https://github.com/soedinglab/MMseqs2.

## Post-Processing


