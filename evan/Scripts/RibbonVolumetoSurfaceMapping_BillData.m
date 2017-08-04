surfoutputdir = '/data/cn5/selfRegulation/V4Process_nosmooth/surf_timecourses/';
fcdatadir = '/data/cn5/selfRegulation/V4Process_nosmooth/final_output/';
surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/';
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
[datapaths subjects ign1 ign2 ign3] = textread(['/data/cn5/selfRegulation/V4Process_nosmooth/FinalDatalist_108.txt'],'%s%s%s%f%f');

HEMS = {'L','R'};

smooth = '2.55';

%subjects = subjects(61:end);

%% Run goodvoxels and map to surface

for s = 1:length(subjects)
%     mkdir(['/data/cn5/selfRegulation/V4Process_nosmooth/final_output/' subjects{s} '/goodvoxels/'])
%     concfilename = ['/data/cn5/selfRegulation/V4Process_nosmooth/final_output/' subjects{s} '/goodvoxels/' subjects{s} '_restruns.conc'];
%     fid = fopen(concfilename,'at'); %open the output file for writing
%     fprintf(fid,'%s\n\r\','number_of_files: 2');
%     fclose(fid);
%     dlmwrite(concfilename,' ','-append');
%     for i = 1:2
%         texttowrite = ['file:' datapaths{s} '/' subjects{s} '/bold' num2str(i) '/' subjects{s} '_b' num2str(i) '_xr3d_333.4dfp.img'];
%         dlmwrite(concfilename,texttowrite,'-append','delimiter','');
%     end
%     
%     disp(['Subject ' num2str(s) ':Run goodvoxels and map to surface'])
%     
%     system(['csh /data/cn4/evan/Scripts/RibbonVolumetoSurfaceMapping_BillData2.csh ' subjects{s} ' ' concfilename])
    
    subfunc = [fcdatadir subjects{s} '/' subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt'];
    system(['niftigz_4dfp -n ' subfunc '.4dfp.img ' subfunc])
    
    submask = ['/data/cn5/selfRegulation/V4Process_nosmooth/final_output/' subjects{s} '/goodvoxels/' subjects{s} '_goodvoxels.nii.gz'];
    
    for hem = 1:2
        
        midsurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.midthickness.native.surf.gii'];
        midsurf_LR32k = [surfdir '/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        whitesurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.pial.native.surf.gii'];
        nativedefsphere = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        outsphere = [surfdir '/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        
        surfname = [subjects{s} '_' HEMS{hem} '_time_333'];
        disp('Map volume to surface')
        system([workbenchdir '/wb_command -volume-to-surface-mapping ' subfunc '.nii.gz ' midsurf ' ' surfoutputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
        disp('Dilate surface timecourse')
        system([workbenchdir '/wb_command -metric-dilate ' surfoutputdir '/' surfname '.func.gii ' midsurf ' 10 ' surfoutputdir '/' surfname '_dil10.func.gii'])
        
        disp('Deform timecourse to 32k fs_LR')
        cd([surfdir '/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/'])
        system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' surfoutputdir '/' surfname '_dil10.func.gii ' surfoutputdir '/' surfname '_dil10_32k_fsLR.func.gii']);
        %system([workbenchdir '/wb_command -metric-resample ' surfoutputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' surfoutputdir '/' surfname '_dil10_32k_fsLR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
        
        disp('Smooth surface timecourse')
        system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' surfoutputdir '/' surfname '_dil10_32k_fsLR.func.gii ' smooth ' ' surfoutputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii'])
        
        system(['rm ' surfoutputdir '/' surfname '.func.gii'])
        system(['rm ' surfoutputdir '/' surfname '_dil10.func.gii'])
        system(['rm ' surfoutputdir '/' surfname '_dil10_32k_fsLR.func.gii'])
    end
    delete([subfunc '.nii.gz'])
end

%% Create smallwall
% disp(['Creating smallwall'])
% for hem = 1:length(HEMS)
%     
%     for s = 1:length(subjects)
%         disp([HEMS{hem} ' hem, subject ' num2str(s)])
%         
%         BOLDname = [surfoutputdir '/' subjects{s} '_' HEMS{hem} '_time_333_dil10_32k_fsLR_smooth' smooth];
%         system(['caret_command64 -file-convert -format-convert ASCII ' BOLDname '.func.gii'])
%         system(['awk ''NF > 25'' ' BOLDname '.func.gii > ' BOLDname '_noHEAD.func.gii'])
%         
%         subdata = load([BOLDname '_noHEAD.func.gii']);
%         subdata = subdata(:,2:end);
%         noproj = single(sum(subdata==0,2) == size(subdata,2));
%         allsubdata(:,s) = noproj;
%         
%         system(['rm ' BOLDname '_noHEAD.func.gii']);
%     end
%     
%     mask = any(allsubdata,2);
%     
%     save(gifti(single(mask)),['/data/cn5/selfRegulation/V4Process_nosmooth/' HEMS{hem} '.atlasroi_group_noproj.func.gii'])
%     
%     save(gifti(single(~mask)),['/data/cn5/selfRegulation/V4Process_nosmooth/' HEMS{hem} '.atlasroi_group_proj.func.gii'])
%     
% end

%% Create ciftis

fcdata_subcortdir = '/data/cn5/selfRegulation/V4Process_nosmooth/final_output_wROIsmooth/';

% ciftioutputdir1 = '/data/cn5/selfRegulation/V4Process_nosmooth/cifti_normalwall_timeseries/';
% mkdir(ciftioutputdir1)
% 
% for s = 1:length(subjects)
%     gunzip([fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii.gz'])
%     disp(['Creating cifti for subject ' num2str(s)])
%     system(['wb_command -cifti-create-dense-timeseries ' ciftioutputdir1 subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii -volume ' fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii /data/cn5/selfRegulation/V4Process_nosmooth/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' surfoutputdir '/' subjects{s} '_L_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii -right-metric ' surfoutputdir '/' subjects{s} '_R_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii'])
%     delete([fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii'])
% end
   

ciftioutputdir2 = '/data/cn5/selfRegulation/V4Process_nosmooth/cifti_smallwall_timeseries_120subcortspace/';
%ciftioutputdir2 = '/data/cn5/selfRegulation/V4Process_nosmooth/cifti_normalwall_timeseries/';
mkdir(ciftioutputdir2)


for s = 1:length(subjects)
    disp(['Subject ' num2str(s)])
    for hem = 1:length(HEMS)

        BOLDname = [surfoutputdir '/' subjects{s} '_' HEMS{hem} '_time_333_dil10_32k_fsLR_smooth' smooth];
        system(['caret_command64 -file-convert -format-convert XML_BASE64 ' BOLDname '.func.gii'])
    end
    gunzip([fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii.gz'])
    disp(['Creating cifti for subject ' num2str(s)])
    %system(['wb_command -cifti-create-dense-timeseries ' ciftioutputdir2 subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii -volume ' fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii /data/cn5/selfRegulation/V4Process_nosmooth/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' surfoutputdir '/' subjects{s} '_L_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/cn5/selfRegulation/V4Process_nosmooth/L.atlasroi_group_proj.func.gii -right-metric ' surfoutputdir '/' subjects{s} '_R_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/cn5/selfRegulation/V4Process_nosmooth/R.atlasroi_group_proj.func.gii'])
    %system(['wb_command -cifti-create-dense-timeseries ' ciftioutputdir2 subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii -volume ' fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' surfoutputdir '/' subjects{s} '_L_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/cn5/selfRegulation/V4Process_nosmooth/L.atlasroi_group_proj.func.gii -right-metric ' surfoutputdir '/' subjects{s} '_R_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/cn5/selfRegulation/V4Process_nosmooth/R.atlasroi_group_proj.func.gii'])
    system(['wb_command -cifti-create-dense-timeseries ' ciftioutputdir2 subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii -volume ' fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' surfoutputdir '/' subjects{s} '_L_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/hcp-zfs/shared-nil/laumannt/120_parcellation/L.atlasroi_group_noproj.func.gii -right-metric ' surfoutputdir '/' subjects{s} '_R_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/hcp-zfs/shared-nil/laumannt/120_parcellation/R.atlasroi_group_noproj.func.gii'])
    %system(['wb_command -cifti-create-dense-timeseries ' ciftioutputdir2 subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii -volume ' fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii /data/cn5/selfRegulation/V4Process_nosmooth/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' surfoutputdir '/' subjects{s} '_L_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii -right-metric ' surfoutputdir '/' subjects{s} '_R_time_333_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii'])
    delete([fcdata_subcortdir subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_smooth2.55.nii'])
end
    
