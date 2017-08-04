%% Create cifti timeseries
TR = 2.5;
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
hemnamelow = {'left';'right'};
smooth = 2.55;
processdir = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/';
outputdir = ['/data/cn4/evan/RestingState/FC_Mapping_120_nocrap/cifti_timeseries_goodvoxels'];
[ign subjects ign2 ign3 ign4] = textread([processdir '/NEW_nokids_DATALIST.txt'],'%s%s%s%f%f');
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';

for s = 1:length(subjects)
    
    subfuncdir = [processdir '/FCPROCESS_bandpass_interp_wROIsmooth/'];
    subfuncvol = [subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55'];
    system(['cp ' subfuncdir '/' subfuncvol '.nii ' outputdir])
    %system(['nifti_4dfp -n ' outputdir '/' subfuncvol '.4dfp.img ' outputdir '/' subfuncvol '.nii'])
   
     surffuncdir = [processdir '/FCPROCESS_bandpass_interp_nosmooth/' subjects{s}];
     surffuncvol = [subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt'];
%      system(['cp ' surffuncdir '/' surffuncvol '.4dfp.* ' outputdir])
%      system(['nifti_4dfp -n ' outputdir '/' surffuncvol '.4dfp.img ' outputdir '/' surffuncvol '.nii'])

        disp(['Processing subject #' num2str(s)])
       
        cd(outputdir)
    
        %Create Cifti timeseries
        disp(['Creating CIFTI timeseries for subject number: ' num2str(s)])
        wmparcdir = '/data/cn4/laumannt/subcortical_mask';
        timename_L = [processdir '/FCPROCESS_bandpass_interp_nosmooth/surf_timecourses/' subjects{s} '_L_time_dil10_32k_fs_LR_smooth' num2str(smooth)];
        timename_R = [processdir '/FCPROCESS_bandpass_interp_nosmooth/surf_timecourses/' subjects{s} '_R_time_dil10_32k_fs_LR_smooth' num2str(smooth)];
        system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_L '.func.gii'])
        system(['caret_command64 -file-convert -format-convert XML_BASE64 ' timename_R '.func.gii'])
        system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth' num2str(smooth) '.dtseries.nii -volume ' subfuncvol '.nii ' wmparcdir '/mode_subcortical_label_LR_333.nii -left-metric ' timename_L '.func.gii -roi-left /data/cn4/evan/RestingState/FC_Mapping_120_nocrap/L.atlasroi_erode3_nocrap.32k_fs_LR.shape.gii -right-metric ' timename_R '.func.gii -roi-right /data/cn4/evan/RestingState/FC_Mapping_120_nocrap/R.atlasroi_erode3_nocrap.32k_fs_LR.shape.gii -timestep ' num2str(TR) ' -timestart 0'])


    system(['rm ' outputdir '/' subfuncvol '*'])
    system(['rm ' outputdir '/' surffuncvol '*'])
end