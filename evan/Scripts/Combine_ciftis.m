%%

subcort120 = load_nii('/data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii');
subcort108 = load_nii('../subcortical_mask/mode_subcortical_LR_333.nii');
subcort120.img = subcort120.img .* single(logical(subcort108.img));
save_nii(subcort120,'mode_subcortical_LR_333.nii')
system('wb_command -volume-label-import mode_subcortical_LR_333.nii /data/cn4/laumannt/subcortical_mask/FreeSurferSubcorticalLabelTableLut_nobrainstem_LR.txt mode_subcortical_label_LR_333.nii -discard-others -drop-unused-labels')
system(['wb_command -cifti-create-dense-timeseries template.dtseries.nii -volume ../final_output_wROIsmooth/01apr12vr_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii.gz mode_subcortical_label_LR_333.nii -left-metric ../surf_timecourses/06mar13rs_L_time_333_dil10_32k_fsLR_smooth2.55.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii -right-metric ../surf_timecourses/06mar13rs_R_time_333_dil10_32k_fsLR_smooth2.55.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii'])
system('wb_command -cifti-correlation template.dtseries.nii corr_template.dconn.nii')
system('wb_command -cifti-convert -to-gifti-ext corr_template.dconn.nii corr_template.func.gii')

%%
bufsize = 524288;

headertext = textread('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_normalwall_timeseries/template.func.gii','%s','delimiter','\r','bufsize',bufsize);
coords108 = zeros(0,3);
for line = 1:size(headertext,1)
    if length(headertext{line}) > 17 && strcmp(headertext{line}(1:17),'<VoxelIndicesIJK>')
        coords108(end+1,:) = str2num(headertext{line}(18:end));
    elseif size(str2num(headertext{line}),2) == 3;
        coords108(end+1,:) = str2num(headertext{line});
    end
end
        
headertext = textread('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/SAIS_116_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii','%s','delimiter','\r','bufsize',bufsize);
coords120 = zeros(0,3);
for line = 1:size(headertext,1)
    if length(headertext{line}) > 17 && strcmp(headertext{line}(1:17),'<VoxelIndicesIJK>')
        coords120(end+1,:) = str2num(headertext{line}(18:end));
    elseif size(str2num(headertext{line}),2) == 3;
        coords120(end+1,:) = str2num(headertext{line});
    end
end 

% %coords108_in120 = zeros(size(coords108,1),1);
% coords108_in120 = [];
% for i = 1:size(coords108,1)
%     if any((coords120(:,1)==coords108(i,1)) .* (coords120(:,2)==coords108(i,2)) .* (coords120(:,3)==coords108(i,3)))
%         reorder = find((coords120(:,1)==coords108(i,1)) .* (coords120(:,2)==coords108(i,2)) .* (coords120(:,3)==coords108(i,3)));
%         coords108_in120(end+1) = reorder;
%     end
% end

%coords120_in108 = zeros(size(coords120,1),1);
coords120_in108 = [];
coords108_in120 = [];
for i = 1:size(coords120,1)
    if any((coords108(:,1)==coords120(i,1)) .* (coords108(:,2)==coords120(i,2)) .* (coords108(:,3)==coords120(i,3)))
        coords120_in108(end+1) = i;
        reorder = find((coords108(:,1)==coords120(i,1)) .* (coords108(:,2)==coords120(i,2)) .* (coords108(:,3)==coords120(i,3)));
        coords108_in120(end+1) = reorder;
    end
end

ncortverts = 29696 + 29716;
inds120 = [(1:ncortverts) (coords120_in108 + ncortverts)];
inds108 = [(1:ncortverts) (coords108_in120 + ncortverts)];


%%


data = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii');
data = data.cdata(inds120,inds120);
data108 = cifti_read('../cifti_correlation_normalwall/108_avg_corr_LR.dconn.nii');
data108 = data108(inds108,inds108);

data = (data .* (120/228)) + (data108 .* (108/228));
%save('avg_corr_120_108_LR','data','-v7.3')
cifti_write_wHDR(data,'corr_template.func.gii','120_108_avg_corr_LR','dconn')

clear data108

cov_corr_L = cov(data(:,1:29696));
save('cov_corr_L','cov_corr_L','-v7.3')
clear cov_corr_L
cov_corr_R = cov(data(:,29697 : (29696 + 29716)));
save('cov_corr_R','cov_corr_R','-v7.3')
clear cov_corr_R



