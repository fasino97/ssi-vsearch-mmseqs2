# ssi-vsearch-mmseqs2
This repository contains all of the files necessary to run accuracy, homogeneity, and completeness using SSI-VSEARCH and SSI-MMseqs2. 

# Databases
The two databases used in the accompanying paper are the Ribosmal Database Project (RDP) and the SILVA rRNA database project:

The RDP database: https://www.glbrc.org/data-and-tools/glbrc-data-sets/ribosomal-database-project.

The SILVA database: https://www.arb-silva.de/.

# Pre-processing
Before you can run any incremental clustering algoirthms you must first use the SplitData.py file to sort the data into testing and training data. VSEARCH expects the input file to be a FASTA file, but MMseqs2 requires formatteed databases that can be converted using MMSeqs2 using createdb.py.

# VSEARCH
All of the information regarding VSEARCH can be found here: https://github.com/torognes/vsearch.

USEARCH is used as well in our SSI-VSEARCH implementation. The download and further information can be found here: https://github.com/rcedgar/usearch12.

# MMSeqs2
The download and inforamtion about MMSeqs2 can be found on this Github repository: https://github.com/soedinglab/MMseqs2.



