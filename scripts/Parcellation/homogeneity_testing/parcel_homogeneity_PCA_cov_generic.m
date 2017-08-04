function output = parcel_homogeneity_PCA_cov_generic(parcellation,cov_corr,datainds)
%outputmetric = parcel_homogeneity_PCA_cov(parcellation,cov_corr,datainds)





parcellation = parcellation(datainds,:);



index = 1;
for map = 1:size(parcellation,2)
    parcellation = parcellation(:,map);
    parcelnum = unique(parcellation);
    parcelnum(parcelnum==0) = [];
    for w = 1:length(parcelnum)
        brainindices{index} = find(parcellation==parcelnum(w));
        mapnum(index) = map;
        eigval_per_first{index} = 0;
        index = index+1;
    end
end

output = zeros(nnz(datainds),1);


%% Calculate PCA, determine explained variance first component

for parcel = 1:length(brainindices)
    
    
    
    parcelcovcorr = cov_corr(brainindices{parcel},brainindices{parcel});
    
    
    if size(parcelcovcorr,1) > 2
        
        eigvals_per = PCA_reduction_cov_onecomp(parcelcovcorr);
        eigval_per_first{parcel} = eigvals_per(1);
        
        
    end
    
end

for parcel = 1:length(eigval_per_first)
    output(brainindices{parcel},mapnum(parcel)) = eigval_per_first{parcel};
end

temp = output;
output = zeros(length(datainds),1);
output(logical(datainds)) = temp;

