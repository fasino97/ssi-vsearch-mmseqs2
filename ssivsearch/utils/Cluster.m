classdef Cluster < handle
    %CLUSTER Cluster class
    
    
    properties 
        %Matrix containing label and count of each tax
        %Column 1 is tax, column 2 is count
        taxCounter = cell(1, 2);
        numLabels = 0;
    end

    
    methods
        %Add new sequence given its tax
        function [] = addSequence(obj, tax)
        %ADDSEQUENCE -  Add tax to cluster
        %ADDSEQUENCE(TAX) increments the counter for a given tax or
        %generates a new tax counter if the tax has not been seen before
            
            %Find index with tax
            %The empty first cell throws off index so add 1 to offset
            indices = strcmp(obj.taxCounter(:,1), tax);
            indices = find(indices(:,1));
            
            
            
            
            %is it already in cluster?
            if(~isempty(indices))
                obj.taxCounter{indices,2} = obj.taxCounter{indices,2} + 1;
            else
                % lastIndex = length(obj.taxCounter) +1;
                obj.numLabels = obj.numLabels + 1;
                
               obj.taxCounter{obj.numLabels, 1} = tax;
               obj.taxCounter{obj.numLabels, 2} = 1;
                
            end
        end
        
        %Get current label of cluster
        function label = getClusterLabel(obj)
            %GETCLUSTERLABEL - Get the label of this cluster
            %GETCLUISTERLABEL() returns the most common taxonomy present in
            %this cluster
            [~, loc] = max([obj.taxCounter{:,2}]);
            if (~isempty(loc))
                label = obj.taxCounter{loc,1};
            else
                label = 'No Label';
            end
        end
        
    end
    
end

