function [ result, binlookup ] = makeTIclasses(grid,numBins)
% this function discretises the distribution of topographic index values
% over a catchment into a user-defined number of classes
% INPUTS:
% - grid: matrix of size equal to the number of grid cells in the domain,
%         containing the topographic index map for the catchment 
%         (i.e. topographic index value for each grid cell)
% - numBins: scalar, equal to the number of classes chosen
% OUTPUTS:
% - result: matrix of size (# of classes,2) containing the TI value for
%           each class in the first column, and the portion of catchment 
%           associated with each class in the second column
% - binlookup: matrix of size equal to the input 'grid', containing, for
%              each grid cell in the catchment, the TI class it belongs to

%%

% reshape grid as an array
array = reshape(grid, size(grid,1)*size(grid,2),1); 

% find which grid cells in array are nans (i.e. outside of catchment area)
nanlookup = isnan(array); %  contains 1 where nan, and 0 otherwise

% remove nans from array
array(nanlookup)=[];

% define TI classes based on the TI values and the chosen number of classes
binEdges = linspace(min(array),max(array), numBins);

% count number of grid cells in array that fall in each class
[h,whichBin] = histc(array, binEdges);

% calculate freq of each class
counts = h / sum(h);

% save TI classes in the first column and freq of each class in the second
% column of a variable 'result'

% column 1 (classes)
tmp = fliplr(binEdges);
result = zeros(numBins,2);
result(:,1) = tmp(1:numBins); 

% column 2 (freq)
tmp = flipud(counts);
result(:,2) = tmp;            

% create variable binlookup
% (as TI classes were defined using the variable array, after having
% removed nans from it, the nans need to be put back in place in order to
% reshape back to a grid -> use of the stored nanlookup)

if(sum(nanlookup) > 0) 
    whichBinWithNan = zeros(length(nanlookup),1);
    index = 1;
    for i = 1:length(nanlookup)
        if(nanlookup(i) == 0) % if grid cell is not nan 
            whichBinWithNan(i) = (numBins - whichBin(index))+1;
            index = index + 1;
        else                  % if grid cell is nan 
            whichBinWithNan(i) = NaN;
        end
    end
else                   
    whichBinWithNan = whichBin;
end

binlookup = reshape(whichBinWithNan, size(grid,1), size(grid,2));

end

