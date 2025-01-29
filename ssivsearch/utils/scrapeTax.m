function [taxonomy] = scrapeTax(taxSequence, depth)
%SCRAPETAX - Find tax at particular depth
%[TAXONOMY] = SCRAPETAX(TAXSEQUENCE, DEPTH) returns the taxonomy at a
% particular taxonomic rank given the full mothur format sequence
%
%   e.g. scrapeTax(Domain;Phylum;Class;Order;Family;Genus, 3)
%           returns Class
%
% taxSequence must not contain the origin Root; prefix and must be between
% 1 and 6, otherwise you will get incorrect answer.
%

taxSequence = split(taxSequence, ';');
if(length(taxSequence) >= depth )
    taxonomy = char(taxSequence(depth));
else
    taxonomy = 'Too Deep';
end

