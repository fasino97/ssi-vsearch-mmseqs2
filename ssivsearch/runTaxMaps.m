%Use the function creatMap(fastaPath,percent) to create a map for data
%taxLookUpMap is for all of the data vs hiddenLabelMap is for testing data
%This step must be done to run the code, but should not affect the clustering

clc;
clear all;
taxLookupMap = createMap('RDP/RDPFull.fasta',1);
hiddenLabelMap = createMap('RDP/RDPtraining.fasta',1);
save('RDPexample.mat', 'taxLookupMap', 'hiddenLabelMap');
