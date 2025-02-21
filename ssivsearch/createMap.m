%createMap function that creates a taxonomic mapping of the sequences to their IDs.
%Percent uses only a percentage of the data in the map, most users will want to map all of the data.
%Use the file runTaxMaps to create a taxLookupMap and a hiddenLabelMap

function [taxMap] = createMap(fastaFile, percent) %enter percent as decimal
    fasta = fastaread(fastaFile);
    createCell = struct2cell(fasta);
    numLines= size(createCell,2);
    
    if(percent == 1)   
        dataList = createCell(1, :);
        tax = cell(1, numLines);
        code = cell(1,numLines);
        for i=1:numLines
            temp = regexp(dataList{i}, 'Root;', 'split');
            code{i} = temp{1}(1:end-1);
            tax{i} = temp{2};
        end
        taxMap = containers.Map(code, tax);
    else
        if(0<percent)&&(percent<1)
        y=numLines*percent;
        ncol = round(y,0);
        x = randperm(size(createCell,2),ncol);
        createCellRand = createCell(:,x);
        numLinesRand=size(createCellRand,2);
        dataListRand = createCellRand(1, :);
        taxRand = cell(1,numLinesRand);
        codeRand = cell(1,numLinesRand);
        for i=1:numLinesRand
            temp = regexp(dataListRand{i}, 'Root;', 'split');
            codeRand{i} = temp{1}(1:end-1);
            taxRand{i} = temp{2};
        end
        taxMap = containers.Map(codeRand, taxRand);
        else
            fprintf('Enter a percent value between 0 and 1\n');
        end  
    end
end
