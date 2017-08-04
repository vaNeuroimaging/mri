cohortfile = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmaskfile = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt'; 

[subjects surfdatafile] = textread(cohortfile,'%s %s');
 [subjects tmasks] = textread(tmaskfile,'%s %s');
 
 workbenchdir = 'env nice -n 4 /data/cn4/laumannt/workbench/bin_linux64/';
 HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
hemname_low = {'left';'right'};
 

for s = 1:length(subjects)
    
    %system(['caret_command64 -file-convert -format-convert XML_BASE64 /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_L_time_dil10_32k_fs_LR_smooth2.55.func.gii /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_R_time_dil10_32k_fs_LR_smooth2.55.func.gii'])
    %system(['caret_command64 -file-convert -format-convert XML_BASE64 /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_2_L_time_dil10_32k_fs_LR_smooth2.55.func.gii /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_2_R_time_dil10_32k_fs_LR_smooth2.55.func.gii'])
    %system(['wb_command -cifti-create-dense-timeseries /data/cn4/evan/RestingState/Ind_variability/Subjects/cifti_timeseries_smallwall/' subjects{s} '_Day1.dtseries.nii -volume /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/ROI_smooth/' subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_wROI255.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_L_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-left /data/hcp-bluearc/home/laumannt/120_parcellation/L.atlasroi_group_noproj.func.gii -right-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY1/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_R_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-right /data/hcp-bluearc/home/laumannt/120_parcellation/R.atlasroi_group_noproj.func.gii -timestep 2.5 -timestart 0'])
    %system(['wb_command -cifti-create-dense-timeseries /data/cn4/evan/RestingState/Ind_variability/Subjects/cifti_timeseries_smallwall/' subjects{s} '_Day2.dtseries.nii -volume /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/ROI_smooth/' subjects{s} '_2_333_zmdt_resid_ntrpl_bpss_zmdt_wROI255.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_2_L_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-left /data/hcp-bluearc/home/laumannt/120_parcellation/L.atlasroi_group_noproj.func.gii -right-metric /data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/DAY2/FCPROCESS_bandpass_nosmooth/surf_timecourses/' subjects{s} '_2_R_time_dil10_32k_fs_LR_smooth2.55.func.gii -roi-right /data/hcp-bluearc/home/laumannt/120_parcellation/R.atlasroi_group_noproj.func.gii -timestep 2.5 -timestart 0'])
    %system(['wb_command -cifti-merge /data/cn4/evan/RestingState/Ind_variability/Subjects/cifti_timeseries_smallwall/' subjects{s} '.dtseries.nii -cifti /data/cn4/evan/RestingState/Ind_variability/Subjects/cifti_timeseries_smallwall/' subjects{s} '_Day1.dtseries.nii -cifti /data/cn4/evan/RestingState/Ind_variability/Subjects/cifti_timeseries_smallwall/' subjects{s} '_Day2.dtseries.nii'])
 
    outputdir = ['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subjects{s} '/'];
    mkdir(outputdir);
    cd(outputdir)

 for hem = 1:2
      %if ~exist([outputdir '/avgcorr_' HEMS{hem} '.dconn.nii'])
         system([workbenchdir '/wb_command -cifti-correlation /data/cn4/evan/RestingState/Ind_variability/Subjects/cifti_timeseries_smallwall/' subjects{s} '.dtseries.nii ' outputdir '/avgcorr_smallwall_' HEMS{hem} '.dconn.nii -roi-override -' hemname_low{hem} '-roi /data/hcp-bluearc/home/laumannt/120_parcellation/' HEMS{hem} '.atlasroi_group_noproj.func.gii -weights ' tmasks{s} ' -fisher-z'])
         disp('Removing NaNs')
           system([workbenchdir '/wb_command -cifti-math "a" avgcorr_smallwall_' HEMS{hem} '_temp.dconn.nii -fixnan 0 -var a avgcorr_smallwall_' HEMS{hem} '.dconn.nii'])
           system(['mv avgcorr_smallwall_' HEMS{hem} '_temp.dconn.nii avgcorr_smallwall_' HEMS{hem} '.dconn.nii'])
           delete(['avgcorr_smallwall_' HEMS{hem} '_temp.dconn.nii'])
         surface_parcellation_cifti_smallwall_wateredge_remover([outputdir '/avgcorr_smallwall_' HEMS{hem} '.dconn.nii'],outputdir,hemname{hem},'yes')
      %end
     
     %surface_parcellation_cifti_smallwall_wateredge_remover([outputdir '/avgcorr_smallwall_' HEMS{hem} '.dconn.nii'],outputdir,hemname{hem},'yes')
     
 end
end
     
     