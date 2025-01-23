classdef IncrementalVsearch < handle
    %INCREMENTALVSEARCH Wrapper around VSearch to incrementalize
    %   Use addBatch with a fasta file to add a collection of sequences.
    %   This assumes the same format as the datasets from mothur. The
    %   header should contain a unique sequence ID and the taxonomy at
    %   each depth separated with a semicolon.
    %
    %   Header Example:
    %       >EU706350_S001045350	Root;Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Nocardioidaceae;Actinopolymorpha;
    %
    %
    %
    %   Use calcPerformance to predict the labels of the unlabeled sequences
    %   and return the number of correct and incorrect predictions.
    %   A cluster that does not have a label is assumed incorrect
    %   but a sequence without a tax with enough depth is not counted
    
    properties
        threshold
        temp
        depth
        taxMap
        hiddenMap
        clusters
        vsearchPath
        usearchPath
        unlabeledSeq
        clusterMembership
        verbose
    end
    
    
    methods
        function [obj] = IncrementalVsearch(threshold, depth, lookupMaps, vsearchPath, usearchPath)
            %INCREMENTALVSEARCH Constructor for IncrementalVsearch
            %INCREMENTALVSEARCH(THRESHOLD, DEPTH, LOOKUPMAPS, VSEARCHPATH,
            %USEARCHPATH) initializes IncrementalVsearch with a temporary
            %folder, a similarity threshold, and a taxonomical depth from
            %1-6.
            %LookupMaps is the name of a .mat file in /utils/ containing
            %a hiddenLabelMap and a taxLookupMap.
            %
            %
            %
            %taxLookupMap should be a containers.Map mapping the sequence
            %ID to the taxonomical depths separated by a semicolon and
            %starting with for all unlabeled and labeled sequences
            %
            % e.g. DQ343153_S000640727 -> Bacteria;Actinobacteria;Actinobacteria;Actinomycetales;Ruaniaceae;Ruania;
            %
            %hiddenLabelMap should be the same as taxLookupMap but only
            %contain labeled sequences
            %vsearchpath and usearchPath are paths to usearch and vsearch
            %and can be generated using the function /utils/getFilePaths.m
            %
            %   Example:
            %       IncrementalVSEARCH(75, 6, 'RDPTaxMaps', ../vsearch.exe, ../usearch.exe)
            %
            
            addpath('utils/');
            addpath(genpath('taxMaps/'));
            obj.threshold = threshold/100;
            obj.depth = depth;
            obj.vsearchPath = vsearchPath;
            obj.usearchPath = usearchPath;
            obj.clusters = Cluster.empty;
            load([lookupMaps, '.mat'] );
            obj.generateID();
            obj.taxMap = taxLookupMap;
            obj.hiddenMap = hiddenLabelMap;
            obj.unlabeledSeq = {};
            obj.clusterMembership = containers.Map;
        end
        
        function [] = generateID(obj)
            %GENERATEID - Generates random temp folder
            %GENERATEID makes a randomly named file using truerand, a script that sends a webrequest to
            %RANDOM.ORG to generate a truely random number based upon
            %atmospheric noise.  MATLAB's rand function I don't believe is
            %acceptable for this application because each IncrementalVSEARCH instance
            %needs its own temp file and there were file collision when running
            %multiple IncrementalVSEARCH instance at the same
            %time on the HPC.  Temp file is of format:
            %'[threshold]_[depth]_[truly random numbers]'
            
            random_nums = truerand(1, 10, 10, 99);
            random_nums = num2str(random_nums);
            n=size(random_nums,1);
            random_nums=random_nums(random_nums~=' ');
            random_nums=reshape(random_nums,n,[]);
            obj.temp = [pwd, '/temp/', num2str(75), '_' , num2str(0), '_' ,random_nums];
            mkdir(obj.temp); %Create the temp directory
            %disp(['Creating temporary folder: ', obj.temp])
        end
        
        
        function [] = addBatch(obj, fastaPath, sort_option)
            %ADDBATCH - Add a fasta file to IncrementalVSEARCH
            %ADDBATCH(FASTAPATH, SORT_OPTION) is the heart of
            %IncrementalVSEARCH It adds a batch of sequences given the path to
            %a fasta file.  It uses the seeds from previous runs stored in the
            %temp folder.  Labels found in hiddenMap are added to clusters,
            %otherwise they are added to the collection of unlabeled sequences
            %sort_option - Use 0 for no sorting, 1 for VSEARCH sorting, 2 for
            %manual sorting
            
            TEMP = obj.temp;
            
            TEMP_FASTA = [TEMP, '/temp.fasta'];
            TEMP_UC = [TEMP, '/temp.uc'];
            TEMP_CENTROIDS = [TEMP, '/temp.centroids'];
            TEMP_OUTCDHIT = [TEMP, '/members.cd-hit'];
            
            if(sort_option ~= 1)
                VSEARCH_RUN = '--usersort --cluster_smallmem ';
            else
                VSEARCH_RUN = '--cluster_fast ';
            end
            
            
            VSEARCH_ARGS = [' --threads 1 ', VSEARCH_RUN, TEMP_FASTA, ' --id ', num2str(obj.threshold), ' --uc ', TEMP_UC, ' --centroids ', TEMP_CENTROIDS];
            %VSEARCH_ARGS = [' --cluster_smallmem ', TEMP_FASTA, ' --id ', num2str(obj.threshold), ' --uc ', TEMP_UC, ' --centroids ', TEMP_CENTROIDS];
            
            
            %Read FASTA
            %disp('Reading fasta file...');
            fasta = fastaread(fastaPath);
            
            %Format FASTA to cell arrray
            fasta_data = cell(length(fasta), 3);
            seqLength = zeros(length(fasta), 1);
            for i=1:length(fasta)
                header = fasta(i).Header;
                taxLoc = strfind(header, 'Root');
                
                tax = header(taxLoc+5:end);
                seqID = header(1:taxLoc-2); %Make sure not to include the R and the tab
                
                fasta_data{i, 1} = tax;
                fasta_data{i, 2} = fasta(i).Sequence;
                fasta_data{i, 3} = seqID;
                seqLength(i) = length(fasta(i).Sequence);
            end
            
            
            %Sort sequences by length
            if(sort_option == 2)
                [~, sortedSeqIndex] = sort(seqLength, 'descend');
                fasta_unsorted_data = fasta_data;
                fasta_data = cell(length(fasta), 3);
                for i = 1:length(fasta_unsorted_data)
                    fasta_data{i, 1} = fasta_unsorted_data{sortedSeqIndex(i), 1};
                    fasta_data{i, 2} = fasta_unsorted_data{sortedSeqIndex(i), 2};
                    fasta_data{i, 3} = fasta_unsorted_data{sortedSeqIndex(i), 3};
                end
            end
            
            
            
            %Read previous runs seeds
            if exist(TEMP_CENTROIDS, 'file') == 2     %Do we have previous seeds saved?
                %disp('Reading previous seeds...')
                previousSeeds = fastaread(TEMP_CENTROIDS);
                currSeeds = cell(length(previousSeeds), 3);
                for i = 1:length(previousSeeds) %Turn seeds to same format as other data
                    % Check if the key exists in the containers.Map object
                    if isKey(obj.taxMap, previousSeeds(i).Header)
                        currSeeds{i, 1} = obj.taxMap(previousSeeds(i).Header);
                        currSeeds{i, 2} = previousSeeds(i).Sequence;
                        currSeeds{i, 3} = previousSeeds(i).Header;
                    else
                    error('Key %s not found in container map.', previousSeeds(i).Header);
                    end
                end
            else
                %disp('No previous seeds found')
                currSeeds = cell(0);
            end
            
            %Generate new fasta with previous seeds and current sequences
            seqIn = vertcat(currSeeds(:,:), fasta_data);
            GenerateFASTA(TEMP_FASTA, seqIn);
            
            %Run VSearch!
            disp(['Running VSEARCH with threshold ', num2str(obj.threshold), ' on ', fastaPath  ,'...'])
            [~,~] = system([obj.vsearchPath, VSEARCH_ARGS]);
            disp('Finished VSEARCH!')
            
            %Convert from UC file type to .clstr
            [~,~] = system([obj.usearchPath, ' --uc2clstr ', TEMP_UC, ' --output ', TEMP_OUTCDHIT]);
            
            
            %Get VSEARCH results
            newBatch = readCDHit(TEMP_OUTCDHIT);
            
            %disp('Reading clusters...')
            %Add sequences to clusters
            for clusterNum = 1:length(newBatch) %For each cluster
                
                if (clusterNum > length(obj.clusters))
                    obj.clusters(clusterNum) = Cluster;
                end
                
                for seqNum = 1:length(newBatch{clusterNum}) %For each sequence
                    %Skip if this sequence was a previous seed
                    if(~isempty(currSeeds))
                        if(ismember(newBatch{clusterNum}{seqNum}, currSeeds(:,3)))
                            continue;
                        end
                    end
                    
                    %Record the cluster number
                    obj.clusterMembership(newBatch{clusterNum}{seqNum}) = clusterNum;
                    
                    if(isKey(obj.hiddenMap,newBatch{clusterNum}{seqNum})) %Is the taxonomy for this sequence known?
                        
                        %Add to cluster
                        thisSequence = obj.hiddenMap(newBatch{clusterNum}{seqNum});
                        thisSequence = scrapeTax(thisSequence, obj.depth);
                        if(~isequal(thisSequence, 'Too Deep'))
                            obj.clusters(clusterNum).addSequence(thisSequence);
                        end
                    else
                        %add tax to unknown collection for later prediction
                        if (clusterNum > length(obj.unlabeledSeq))
                            obj.unlabeledSeq{clusterNum} = {};
                        end
                        
                        obj.unlabeledSeq{clusterNum}{end+1} = newBatch{clusterNum}{seqNum};
                    end
                end
            end
            
            
            disp('Done!')
            %disp(' ')
        end
        
        
        function [numCorrect, numIncorrect]  = calcPerformance(obj)
            %CALCPERFORMANCE - Calculate the classification performance on unlabeled
            %sequences.  A cluster that does not have a label or a sequence without a tax with enough depth is not counted.
            numCorrect = 0;
            numIncorrect = 0;
            for clusterNum = 1:length(obj.unlabeledSeq)
                thisClusterLabel = obj.clusters(clusterNum).getClusterLabel(); %The mode of the cluster tax
                if(~strcmp(thisClusterLabel, 'No Label'))
                    numIncorrect = numIncorrect + 1; %was empty (which would mean no labels aren't counted)
                    for seqNum = 1:length(obj.unlabeledSeq{clusterNum})
                        if(~isempty(obj.unlabeledSeq{clusterNum}{seqNum}))
                            thisSeqTax = obj.taxMap(obj.unlabeledSeq{clusterNum}{seqNum}); %Tax of unlabeled sequence
                            thisSeqTax = scrapeTax(thisSeqTax, obj.depth);
                            if(strcmp(thisClusterLabel, thisSeqTax)) %Is the prediction correct?
                                numCorrect = numCorrect + 1;
                            else
                                if(~strcmp(thisSeqTax, 'Too Deep')) %Only punish if we have the depth
                                    numIncorrect = numIncorrect + 1;
                                end
                            end
                            
                        end
                    end   
                end
            end
            
            obj.unlabeledSeq = {}; %Remove all the unlabeled sequences
        end
    end
    
end

