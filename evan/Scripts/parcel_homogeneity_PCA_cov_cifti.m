function output = parcel_homogeneity_PCA_cov_cifti(parcels,cov_corr)



output = zeros(size(parcels));

parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];



%% Calculate PCA, determine explained variance first component
for parcelnum = 1:length(parcelIDs)
    inds = (parcels==parcelIDs(parcelnum));
        parcelcovcorr = cov_corr(inds,inds);
        
        
        if size(parcelcovcorr,1) > 2
            
            eigvals_per = PCA_reduction_cov_onecomp(parcelcovcorr);
            output(inds) = eigvals_per(1);
            
        end
end

    
