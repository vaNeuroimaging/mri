function overlap = WatershedAssign(watershed1name, watershed2name)
% watershed1name = C1;
% watershed2name = C2;

if ischar(watershed1name)
    watershed1 = gifti(watershed1name); watershed1 = watershed1.cdata;
else
    watershed1 = watershed1name;
end

if ischar(watershed2name)
    watershed2 = gifti(watershed2name); watershed2 = watershed2.cdata;
else
    watershed2 = watershed2name;
end

watershed1parcelIDs = unique(watershed1); watershed1parcelIDs(watershed1parcelIDs==0) = [];
watershed2parcelIDs = unique(watershed2); watershed2parcelIDs(watershed2parcelIDs==0) = [];

costmat = ones(length(watershed1parcelIDs),length(watershed2parcelIDs));

for i = 1:length(watershed1parcelIDs)
    for j = 1:length(watershed2parcelIDs)
        
        numoverlap = nnz((watershed1==watershed1parcelIDs(i)) .* (watershed2==watershed2parcelIDs(j)));
        
        costmat(i,j) = mean([(numoverlap/nnz(watershed1==watershed1parcelIDs(i))) , (numoverlap/nnz(watershed2==watershed2parcelIDs(j)))]);
        
    end
end

[assignment, cost] = munkres((1 - costmat));

overlap = zeros(size(watershed1));

for i = 1:length(watershed1parcelIDs)
    
    if assignment(i) > 0 && costmat(i,assignment(i)) > 0;
        
        indices = logical((watershed1==watershed1parcelIDs(i)) .* (watershed2==watershed2parcelIDs(assignment(i))));
        overlap(indices) = costmat(i,assignment(i));
    end
end
        