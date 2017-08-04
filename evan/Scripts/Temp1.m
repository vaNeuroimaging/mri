avgcrosscorr(isnan(avgcrosscorr)) = 0;
covmat = cov(avgcrosscorr);

for parcelcounti = 1:length(allrandparcels)
    for iter = 1:100
        string{iter} = ['Threshold ' num2str(parcelcounti) ', Iteration ' num2str(iter)];
        if iter==1; fprintf('%s',string{iter}); else fprintf([repmat('\b',1,length(string{iter-1})) '%s'],string{iter}); end
    
        
        label = allrandparcels{parcelcounti}(:,iter);
        
        PCAvals = parcel_homogeneity_PCA_vector(label,covmat,'L');
   
           watersheds = unique(label);
   watersheds(watersheds==0) = [];
   newlabel = zeros(size(label));
   for watershednum = 1:length(watersheds)
       watershed = watersheds(watershednum);
       newlabel(label==watershed) = watershednum;
       eigvals_per_first{parcelcounti}(watershednum,iter) = mean(PCAvals(newlabel==watershednum));
   end
    
   multilabel(:,iter) = newlabel;

end

disp(' ')



newallrandparcels{parcelcounti} = multilabel;

end