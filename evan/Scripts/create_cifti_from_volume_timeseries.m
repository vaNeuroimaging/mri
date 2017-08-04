function create_cifti_from_volume_timeseries(volume,outputfile)

if strcmp(volume(end-8:end),'.4dfp.img')
    system(['nifti_4dfp -n ' volume ' ' volume(1:end-9) '.nii'])
    volume = [volume(1:end-9) '.nii'];
end

map_vol_to_surface(volume)

if ~exist(outputfile)
    outputfile = [volume(1:end-4) '_eroded.dtseries.nii'];
end

system(['wb_command -cifti-create-dense-timeseries ' outputfile ' -volume ' volume ' /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' volume(1:end-4) '_L_32k_fsLR.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii -right-metric ' volume(1:end-4) '_R_32k_fsLR.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii'])