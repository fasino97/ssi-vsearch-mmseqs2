function [] = GenerateFASTA(filename, input)
%GENERATEFASTA - GENERATE A FASTA FILE IN PROPER MOTHUR
%FORMAT
%GENERATE(FILENAME, INPUT) creates and writes to a fasta file that is in
%the proper mothur format that can be read with IncrementalVSEARCH.  The
%input is a cell array.

%   Format of 'input' :
%   col 1: Sequence Labels -> Root;Bacteria;etc..
%   col 2: Actual RNA Sequence -> TCGACGT...
%   col 3: Sequence ID 


% currentLine = 1;
fileID = fopen(filename,'w'); % get a write handle

for n = 1:length(input)
    fprintf(fileID, '>%s\tRoot;%s\r\n', input{n, 3}, input{n,1}); % Header line
    fprintf(fileID, '%s\r\n', input{n, 2}); % Sequence line
end

fclose(fileID); % close our write handle

end