tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');
%subjects = subjects(1:120);
%tmasks = tmasks(1:120);

ngrayords = 66697;

dtseries_template = [];
dconn_template = '/data/hcp-zfs/shared-nil/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii';

datapresent = zeros(ngrayords,1);
corr = zeros(ngrayords,ngrayords,2);
for s = 1:length(subjects)
    subname = subjects{s};
    disp(s)
    datafile = ['/data/hcp-zfs/shared-nil/laumannt/120_parcellation/FCPROCESS_ZEROPAD/FCPROCESS_SCRUBBED/cifti_timeseries_normalwall/' subname '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = cifti_read(datafile);
    tmask = load(tmasks{s});
    data = data(:,logical(tmask));
    datapresent = datapresent + ~all(data==0,2);
    corr(:,:,2) = FisherTransform(paircorr_mod(data'));
    corr(:,:,1) = nansum(corr,3);
end
corr = corr(:,:,1);

cifti_write_wHDR(datapresent / length(subjects),dtseries_template,'120_datapresent_LR.dtseries.nii')

corr = corr / length(subjects);
cifti_write_wHDR(corr,dconn_template,'120_avg_corr_LR.dconn.nii','dconn')
    