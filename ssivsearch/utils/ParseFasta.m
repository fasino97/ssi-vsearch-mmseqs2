function [fasta_data] = ParseFasta(fasta_filename)

fasta = fastaread(fasta_filename);


fasta_data = cell(length(fasta), 3);
for i = 1:length(fasta)
    header = textscan(fasta(i).Header, '%s');
    seqID = header{1}{1};
    tax = header{1}{2}(6:end-1);
    
    
    fasta_data{i, 1} = tax;
    fasta_data{i, 2} = fasta(i).Sequence;
    fasta_data{i, 3} = seqID;
    
    
end


end