[subjects tmasks] = textread(['/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_108.txt'],'%s%s');
datafolder = '/data/cn5/selfRegulation/V4Process_nosmooth/cifti_smallwall_timeseries/';
outfolder = '/data/cn5/selfRegulation/V4Process_nosmooth/cifti_correlation_smallwall/';
dconntemplate = '/data/cn5/selfRegulation/V4Process_nosmooth/cifti_correlation_smallwall/template.func.gii';

Lmask = gifti('/data/cn5/selfRegulation/V4Process_nosmooth/L.atlasroi_group_proj.func.gii'); Lmask = Lmask.cdata;
Rmask = gifti('/data/cn5/selfRegulation/V4Process_nosmooth/R.atlasroi_group_proj.func.gii'); Rmask = Rmask.cdata;

Lnwmask = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'); Lnwmask = ~Lnwmask.cdata;
L_maskinmask = Lnwmask(logical(Lmask));
Rnwmask = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'); Rnwmask = ~Rnwmask.cdata;
R_maskinmask = Rnwmask(logical(Rmask));



for s = 1:length(subjects)
    disp(['Subject ' num2str(s)])
    tmask = load(tmasks{s});
    timeseries = cifti_read([datafolder subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii']);
    timeseries = timeseries(:,logical(tmask))';
    timeseries(isnan(timeseries)) = 0;
    
    correl = FisherTransform(paircorr_mod(timeseries)) ./ length(subjects);
    correl(isnan(correl)) = 0;
    if s==1
        avgcorrel = correl;
    else
        avgcorrel = avgcorrel + correl;
    end
    clear correl
end
cd(outfolder)
cifti_write_wHDR(avgcorrel,dconntemplate,['108_avg_corr_LR'],'dconn')  



hemcorrel = avgcorrel(1:nnz(Lmask),:);
cifti_write_wHDR(hemcorrel,'/data/cn5/selfRegulation/V4Process_nosmooth/cifti_correlation_smallwall/templateL.func.gii',['108_avg_corr_L'])

all_maskinmask = [L_maskinmask;R_maskinmask;ones(size(hemcorrel,2) - length(L_maskinmask) - length(R_maskinmask) , 1)];

hemcorrel = hemcorrel(logical(L_maskinmask),logical(all_maskinmask));
cov_corr_L = cov(hemcorrel');
save('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_covariance_normalwall/cov_corr_L_normalwall.mat','cov_corr_L','-v7.3')
clear hemcorrel cov_corr_L



hemcorrel = avgcorrel((nnz(Lmask)+1):(nnz(Lmask)+nnz(Rmask)),:);
cifti_write_wHDR(hemcorrel,'/data/cn5/selfRegulation/V4Process_nosmooth/cifti_correlation_smallwall/templateR.func.gii',['108_avg_corr_R'])

hemcorrel = hemcorrel(logical(R_maskinmask),logical(all_maskinmask));
cov_corr_R = cov(hemcorrel');
save('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_covariance_normalwall/cov_corr_R_normalwall.mat','cov_corr_R','-v7.3')
clear hemcorrel cov_corr_R
