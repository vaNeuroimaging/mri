testsub = 5;

counter = 0;
for s = 1:length(subjects)
    FAfile = ['../subjects/' subjects(s).name '/DTI/FA_near_cortex_LR.dscalar.nii'];
    if exist(FAfile)
        counter = counter+1;
        FA = ft_read_cifti_mod(FAfile);
        all_FA(:,:,counter) = FA.data;
    end
end
%FA_avg = mean(all_FA,3);
%FA.data = FA_avg;
%ft_write_cifti_mod('FA_near_cortex_avgMAV',FA)

tvals = zeros(size(FA.data));
subvec = setdiff([1:size(all_FA,3)],testsub);
for i = 1:numel(FA.data);
    [x,y] = ind2sub(size(FA.data),i);
    [H,P,CI,STATS] = ttest(all_FA(x,y,subvec),all_FA(x,y,testsub));
    
    tvals(i) = STATS.tstat;
    
end

FA.data = tvals;
ft_write_cifti_mod('Ttests_FA_near_cortex_MAV006_v_otherMAV',FA)
    