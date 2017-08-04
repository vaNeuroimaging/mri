%volume is Temp.nii

for hem = 1:2
    
    surfname = ['Temp_' HEMS{hem}];
    
    
    midsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.midthickness.native.surf.gii'];
    midsurf_LR32k = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    whitesurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.white.native.surf.gii'];
    pialsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.pial.native.surf.gii'];
    nativedefsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
    outsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
    
    disp('Map volume to surface')
    system([workbenchdir '/wb_command -volume-to-surface-mapping ' outputdir '/Temp.nii ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
    disp('Dilate surface timecourse')
    system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
    
    disp('Deform timecourse to 32k fs_LR')
    cd([surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/'])
    %system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' outputdir '/' surfname '_dil10.func.gii ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii']);
    system([workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
    
    disp('Smooth surface timecourse')
    system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii ' smooth ' ' outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii'])
    
end

disp('Creating dense timeseries')
    system(['wb_command -cifti-create-dense-timeseries ' outputdir '/' subname{subject} '.dtseries.nii -volume ' outputdir '/Temp.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' outputdir '/Temp_L_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii -right-metric ' outputdir '/Temp_R_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii'])