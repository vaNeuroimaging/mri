subnames = dir('/home/data/subjects/vc*');

for i = 1:length(subnames)
    subject = subnames(i).name;
    data = ft_read_cifti_mod(['/home/data/subjects/' subject '/fs_LR/NativeVol/fsaverage_LR32k/' subject '.LR.thickness.32k_fs_LR.dtseries.nii']);
    if i==1
        all = data;
    else
        all.data(:,i) = data.data;
    end
end
ft_write_cifti_mod('DART_cortical_thickness',all)
all.data = mean(all.data,2);
ft_write_cifti_mod('DART_cortical_thickness_mean',all)