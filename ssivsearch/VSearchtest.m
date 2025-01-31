clear
clc

%File used to run incremental VSEARCH. 
%The example below shows a 97 percent similarity threshold for five training batches and one testing batch.
%phylum97 is a place holder variable that can be named anything.
%(similarity threshold, taxonomic depth, associated map, vsearch.exe, usearch.exe)
%Similarity treshold ranges from 75-97.
%Taxonomic depth does not matter for just running VSEARCH, only for the included accuracy metric.
%Map can be created using createMap function

phylum97=IncrementalVsearch(97,2,'Silva2024map','vsearch.exe','usearch.exe');
phylum97.addBatch('Silva/silvatraining2024batch1.fasta',2);
phylum97.addBatch('Silva/silvatraining2024batch2.fasta',2)
phylum97.addBatch('Silva/silvatraining2024batch3.fasta',2)
phylum97.addBatch('Silva/silvatraining2024batch4.fasta',2)
phylum97.addBatch('Silva/silvatraining2024batch5.fasta',2)
phylum97.addBatch('Silva/silvatesting2024.fasta',2)

