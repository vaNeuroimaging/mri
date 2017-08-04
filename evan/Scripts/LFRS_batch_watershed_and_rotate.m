cohortfile = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmaskfile = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt'; 

[subjects surfdatafile] = textread(cohortfile,'%s %s');
 [subjects tmasks] = textread(tmaskfile,'%s %s');
 
 workbenchdir = 'env nice -n 4 /data/cn4/laumannt/workbench/bin_linux64/';
 HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
hemname_low = {'left';'right'};
 

for s = 1:length(subjects)
    
    
    cd(['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subjects{s} '/'])
    if exist(['avgcorrofcorr_allgrad_L_smooth2.55_wateredge_thresh0.01_avg.func.gii'])
        
        
         system(['wb_command -cifti-create-dense-timeseries ' subjects{s} '_Day1_normalwall.dtseries.nii -volume /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/ROI_smooth/' subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_wROI255.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_L_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii -right-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_R_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii -timestep 2.5 -timestart 0'])
            system(['wb_command -cifti-create-dense-timeseries ' subjects{s} '_Day2_normalwall.dtseries.nii -volume /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/ROI_smooth/' subjects{s} '_2_333_zmdt_resid_ntrpl_bpss_zmdt_wROI255.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_2_L_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii -right-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_2_R_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii -timestep 2.5 -timestart 0'])
            system(['wb_command -cifti-merge ' subjects{s} '_timeseries_normalwall.dtseries.nii -cifti ' subjects{s} '_Day1_normalwall.dtseries.nii -cifti ' subjects{s} '_Day2_normalwall.dtseries.nii'])
            delete([subjects{s} '_Day*_normalwall.dtseries.nii'])
           
        
        
       for hemnum = 1:length(HEMS)
           hem = HEMS{hemnum};
           
           watershed_algorithm_merge_andthresh_perc(['avgcorrofcorr_allgrad_' hem '_smooth2.55_wateredge_thresh0.01_avg.func.gii'],'./',[subjects{s} '_' hem '_'],hem,.45);
           
           system(['wb_command -cifti-correlation ' subjects{s} '_timeseries_normalwall.dtseries.nii ' subjects{s} '_corr_' hem '_normalwall.dtseries.nii -roi-override -' hemname_low{hemnum} '-roi /data/cn4/laumannt/subcortical_mask/' hem '.atlasroi.32k_fs_LR.shape.gii -weights ' tmasks{s} ' -fisher-z'])
           disp('Removing NaNs')
           system([workbenchdir '/wb_command -cifti-math "a" ' subjects{s} '_corr_' hem '_normalwall_temp.dtseries.nii -fixnan 0 -var a ' subjects{s} '_corr_' hem '_normalwall.dtseries.nii'])
           system(['mv ' subjects{s} '_corr_' hem '_normalwall_temp.dtseries.nii ' subjects{s} '_corr_' hem '_normalwall.dtseries.nii'])
           system(['wb_command -cifti-convert -to-gifti-ext ' subjects{s} '_corr_' hem '_normalwall.dtseries.nii ' subjects{s} '_corr_' hem '_normalwall.func.gii'])
    
           delete([subjects{s} '_corr_' hem '_normalwall.dtseries.nii'])
           
           a = gifti([subjects{s} '_corr_' hem '_normalwall.func.gii']); %a = a.cdata(:,(1:nnz(mask))+(strcmp(hem,'R') * ncortLverts));

            cov_corr = cov(a.cdata');
            clear a
            save(['cov_corr_' hem '.mat'],'cov_corr','-v7.3')
            
            generate_rotated_parcels_andPCA6([subjects{s} '_' hem '_watershedmerge_0.45.func.gii'],1000,cov_corr,1,hem,[subjects{s} '_' hem '_'])
            generate_rotated_parcels_andPCA6(['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45.func.gii'],1000,cov_corr,1,hem,['Group_' hem '_'])
%            generate_rotated_parcels_andPCA6(['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh2_watershedmerge_0.45.func.gii'],1000,cov_corr,1,hem,['Group_' hem '_'])
            
       end
    end
end