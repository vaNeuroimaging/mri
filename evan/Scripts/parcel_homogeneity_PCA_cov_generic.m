function outputmetric = parcel_homogeneity_PCA_cov_generic(watershed,cov_corr,datainds)
%outputmetric = parcel_homogeneity_PCA_cov(watershed,cov_corr,iscifti,hem)





        watershed = watershed(datainds,:);
    


index = 1;
for map = 1:size(watershed,2)
    watershed = watershed(:,map);
    waternum = unique(watershed);
    waternum(waternum==0) = [];
    for w = 1:length(waternum)
        brainindices{index} = find(watershed==waternum(w));
        mapnum(index) = map;
        eigval_per_first{index} = 0;
        index = index+1;
    end
end

outputmetric = zeros(nnz(datainds),1);


%% Calculate PCA, determine explained variance first component

for parcel = 1:length(brainindices)
    

    
        watercovcorr = cov_corr(brainindices{parcel},brainindices{parcel});
        
        
        if size(watercovcorr,1) > 2
            
            eigvals_per = PCA_reduction_cov_onecomp(watercovcorr);
            eigval_per_first{parcel} = eigvals_per(1);
            

        end
    
end

for parcel = 1:length(eigval_per_first)
    outputmetric(brainindices{parcel},mapnum(parcel)) = eigval_per_first{parcel};
end
    
temp = outputmetric;
outputmetric = zeros(length(datainds),1);
outputmetric(logical(datainds)) = temp;

