%Reads .cd-hit file and returns clusters
function [clusters] = readCDHit(filename)

input = fastaread(filename);
numClusters = length(input);
clusters = cell(numClusters, 1);

% for each cluster
for n=1:numClusters
    splitCell = split(input(n).Sequence, ',');
    
    for i=2:length(splitCell)
        
        thisString = char(splitCell(i));
        
        if(contains(thisString, 'Root'))
            clusters{n}{i-1} = thisString((strfind(thisString, '>')+1):(strfind(thisString, 'Root')-1)); %Extract ID
        else
            clusters{n}{i-1} = thisString((strfind(thisString, '>')+1):(strfind(thisString, '...')-1)); %Extract ID
        end
        
    end
    
end