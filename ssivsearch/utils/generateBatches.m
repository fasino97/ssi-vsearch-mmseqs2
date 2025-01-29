function [  ] = generateBatches( fastaPath, seqPerBatch )
%GENERATEBATCHES Summary of this function goes here
%   Detailed explanation goes here

fasta = fastaread(fastaPath);
numSeq = length(fasta);


% make an array with {tax, seq}, pre-allocate for speed
unsortedSeq = cell(numSeq, 2); % un-sorted
seq = cell(numSeq, 2); % sorted
seqLength = zeros(numSeq, 1);

for i=1:numSeq
    header = fasta(i).Header;
    taxLoc = strfind(header, 'Root');
    tax = header(taxLoc:end);
    seqID = header(1:taxLoc);
    
    unsortedSeq{i, 1} = seqID;
    unsortedSeq{i, 2} = tax;
    
    
    seqLength(i) = length(fasta(i).Sequence);
end

% sort them from largest -> smallest
% maintain the tax info!!!
[sortedSeqLengths, sortedSeqIndex] = sort(seqLength, 'descend');

% populate our new cell array w/ sorted seq
% toss in the header in column 3, as well.
for n=1:numSeq
    seq{n, 1} = unsortedSeq{sortedSeqIndex(n), 1};
    seq{n, 2} = unsortedSeq{sortedSeqIndex(n), 2};
end


%%%%%%%UNFINISHED%%%%%%

% Shuffle
%seq = shuffleVector(seq);


% slice into multiple batches
numBatches = ceil(numSeq / seqPerBatch);

batch = cell(numBatches, 1);
% col 1: Sequence ID
% col 2: Tax
for n=1:numBatches
    
    indBegin = (n-1)*seqPerBatch + 1;
    indEnd = n*seqPerBatch;
    
    % Bounds checking - our last batch may not be a full seqPerBatch long.
    if (indEnd > numSeq)
        indEnd = numSeq;
    end
    
    batch{n} = seq(indBegin:indEnd, :);
    
end


end

