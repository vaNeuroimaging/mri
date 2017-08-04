function output_indices = matrix_thresholder_faster2_structurespecific_ramefficient(matrix,threshold,structinds)
%
% output_indices = matrix_thresholder_faster2_structurespecific_ramefficient(matrix,threshold,structinds)
%
%EMG modified, 03/2015

d=size(matrix);

%use upper triangle of matrix
matrix = matrix(triu(true(d),1));

%get structural divisions
structs = unique(structinds);

divisioncounter = 0;

%get indices in matrix corresponding to pairings of structures
inds = ones(d,'uint8');
for structnum1 = 1:length(structs)
    for structnum2 = structnum1:length(structs)
        divisioncounter = divisioncounter +1;
        inds(structinds==structs(structnum1),structinds==structs(structnum2)) = divisioncounter;
    end
end
%get upper triangle of indices
inds = inds(triu(true(d),1));

%for each pairing of structures
for division = 1:divisioncounter
    
    %get the indices of this pairing of structures
    smallinds = find(inds==division);
    
    %get the matrix values for this pairing of structures
    smallvals = matrix(smallinds);
    
    %sort those values
    [~, smallorder]=sort(smallvals,'descend');
    clear smallvals
    
    %put the top (kden) values back into the matrix as ones and the rest as zeros
    matrix(smallinds(smallorder)) = [true(ceil(numel(smallorder) * kden_thresh),1); false(floor(numel(smallorder) * (1-kden_thresh)),1)];
    clear smallorder smallinds
    
end
clear inds
matrix = logical(matrix);

%turn the upper triangle thresholded values back into a full matrix
matrixout = false(d); matrixout(triu(true(d),1)) = matrix;
clear matrix

%get the indices of cells that were above the threshold
output_indices = find(matrixout);





