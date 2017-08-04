subjectnames = {'Movie006'};
for s = 1:length(subjectnames)

subjectname = subjectnames{s};

RSFC = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/connectome/RSFC_corr.dconn.nii']);
Movie = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/connectome/Movie_corr.dconn.nii']);

out = zeros(length(RSFC.data),1);
for i = 1:length(RSFC.data)
    out(i) = paircorr_mod(RSFC.data(:,i),Movie.data(:,i));
end

RSFC.data = out;
RSFC.dimord = 'pos_time';
ft_write_cifti_mod([subjectname '_RSFC_v_Movie_correlation'],RSFC)
clear Movie
end