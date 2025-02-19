# ssi-vsearch-mmseqs2
This repository contains all of the files necessary to run accuracy, homogeneity, and completeness using SSI-VSEARCH and SSI-MMseqs2. 

## Databases and Necessary Downloads
The two databases used in the accompanying paper are the Ribosmal Database Project (RDP) and the SILVA rRNA database project:

The RDP database: https://www.glbrc.org/data-and-tools/glbrc-data-sets/ribosomal-database-project.

The SILVA database: https://www.arb-silva.de/.

Vsearch, Usearch and MMSeqs2 are all necessary to make all aspects of this repository work.
All of the information regarding VSEARCH can be found here: https://github.com/torognes/vsearch.

USEARCH is used as well in our SSI-VSEARCH implementation. The download and further information can be found here: https://github.com/rcedgar/usearch12.

The download and inforamtion about MMSeqs2 can be found on this Github repository: https://github.com/soedinglab/MMseqs2.

## Pre-Processing
Before you can run any incremental clustering algoirthms you must first use the SplitData.py file to sort the data into testing and training data. VSEARCH expects the input file to be a FASTA file, but MMseqs2 requires formatteed databases that can be converted using MMSeqs2 using createdb.py.

## SSI-VSEARCH
SSI-VSEARCH is run in Matlab. The main code is run in VSearchtest.m, where you specify the similarity threshold (75-97), taxonomic depth (2-6), associated map, and the vsearch and usearch executable files. Taxonomic depth does not need to be set for running SSI-VSEARCH, as this is only used in the accuracy calculations available when running the IncrementalVSearch.m file. For this improved version of SSI-VSEARCH we use a more specific post-processing code. The associated map can be created by running the createMap function where the pecentage of the dataset being mapped can be set. Most users will want to map the entire training dataset. 

## SSI-MMSeqs2
SSI-MMSeqs2 is run in Python. This code was run mainly on Rowan University's HPC so the code attached (slurm-script) is the code to call MMSeqs2 and capture the time it takes to run. A database must be used as the input by using MMSeqs2's createdb function, and batches are added using MMSeqs2's clusterupdate function.

Example of running SSI-MMSeqs2:
1. ```mmseqs createdb training.fasta trainingdb```
2. ```mmseqs createdb trainingplustesting.fasta trainingplustestingdb```
3. ```mmseqs cluster trainingdb trainingclusters tmp```
4. ```mmseqs clusterupdate trainingdb trainingplustestingdb trainingclusters newcombineddb trainingplustestingclusters tmp```
5. ```mmseqs create tsv newcombineddb newcombineddb trainingplustestingclusters clusters.tsv```
Tags can be included at the end of those commands to specific similarity threshold. Examples used in this experiment include: ```--min-seq-id 0.75``` tp set the similarity threshold to 75%, ```--cluster-mode 0``` to make sure the clustering is the default method.

## Post-Processing
Post-processing for SSI-VSEARCH and SSI-MMSeqs2 can be found in two seperate codes. SSI-VSEARCH's post-processing code requires .cd-hit files. SSI-MMSeqs2's post-processing requires a .tsv file. This is not the file immeditately given when MMSeqs2 clusters, so that output result must be put through MMSeqs2's createtsv function to get the proper file format. Both codes run similarly; inputs for both are a training fasta file, a testing fasta file, and the clustering output file. 

