
PCAresults = {'PCA_eigval_per_first_water_Smoothedge_LR_minima5_watershed.dtseries.nii' 'PCA_eigval_per_first_water_Edge_LR_minima5_watershed.dtseries.nii'};

finalvals = [];
finalsize = [];

for i = 1:length(PCAresults)

    PCAdata{i} = cifti_read(PCAresults{i});

    uniquevals{i} = unique(PCAdata{i});
    uniquevals{i}(uniquevals{i}==0) = [];
    
    for parcel = 1:length(uniquevals{i});
        
        parcelsizes{i}(parcel,:) = length(find(PCAdata{i}==uniquevals{i}(parcel)));
        
    end
    
    finalvals = [finalvals; uniquevals{i}];
    finalsize = [finalsize; parcelsizes{i}];
    
end

[B Bint Resid] = regress(finalvals, [finalsize ones(length(finalsize),1)]);

residcounter = 1;

for i = 1:length(PCAresults)
    
    outputdata{i} = zeros(size(PCAdata{i}));
    
    for parcel = 1:length(uniquevals{i});
        
        parcelindices = find(PCAdata{i}==uniquevals{i}(parcel));
        
        outputdata{i}(parcelindices) = Resid(residcounter) + B(2);
        
        residcounter = residcounter+1;
        
    end
    
    cifti_write(outputdata{i},[PCAresults{i}(1:end-13) '_sizeregressed'])
    
end

