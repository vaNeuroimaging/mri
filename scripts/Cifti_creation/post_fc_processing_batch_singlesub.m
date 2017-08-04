function systemresult = post_fc_processing_batch_singlesub(subject,voldata_tomap,outfolder,fsLRfolder,TR,sequence_name,tmaskfile,preproc_datalist,runs_sessionsfile,medial_masks,smoothnum,systemresult)
%systemresult = post_fc_processing_batch_singlesub(subject,voldata_tomap,TR,sequence_name,tmaskfile,preproc_datalist,post_fc_processing_batch_params_file,systemresult)
%
% Run through post fc-processing, creating cifti timeseries from
% unsmoothed fc-processed functional data.
%

if ~exist('systemresult')
    systemresult = cell(0,2);
end


subcort_mask = '/home/data/scripts/Resources/cifti_masks/subcortical_mask_LR_333_MNI.nii.gz';

medial_mask_L = medial_masks{1,1};
medial_mask_R = medial_masks{1,2};


if size(medial_masks,1) == 1
    do_smallwall = false;
else
    do_smallwall = true;
    sw_medial_mask_L = medial_masks{2,1};
    sw_medial_mask_R = medial_masks{2,2};
end



%-------------------------------------------------------------------------

%fsLRfolder = ['/home/data/subjects/' subject '/fs_LR/MNI/'];

%outfolder = ['/home/data/subjects/' subject '/cifti/'];

warning off
mkdir(outfolder)
surffuncdir = [outfolder '/surf_timecourses/'];
ciftidir = [outfolder '/cifti_timeseries_normalwall/'];
ciftiswdir = [outfolder '/cifti_timeseries_smallwall/'];
goodvoxfolder = [outfolder '/goodvoxels/'];
outputdatalistname = [outfolder '/cifti_datalist.txt'];
outputdatalistname_sw = [outfolder '/cifti_datalist_smallwall.txt'];

workbenchdir = '/usr/local/workbench/bin_rh_linux64/';
HEMS = {'L';'R'};

prevstring = [];

%Make folders
mkdir(surffuncdir);
mkdir(ciftidir);
mkdir(goodvoxfolder);
if do_smallwall
    mkdir(ciftiswdir);
end


%Convert and copy data over
%disp(['Subject ' subject ': copying functional data']);

if ~exist([voldata_tomap '.nii.gz']) && ~exist(voldata_tomap)
    error(['subject data ' voldata_tomap ' does not exist!'])
end
funcvol = [surffuncdir '/' sequence_name '_funcvol'];

%Remove NaNs from the data and copy to new location
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' voldata_tomap ' -nan ' funcvol]);


%Identify high-SNR voxels
%disp(['Subject ' subject ': calculating high-SNR voxels from unprocessed functional data']);

if ~isempty(preproc_datalist)
    
    [preproc_runs,~,~] = textread(preproc_datalist,'%s%s%s');
    [runIDs,sessIDs] = textread(runs_sessionsfile,'%n%n');
    
    tmask = load(tmaskfile);

    for runnum = 1:length(preproc_runs)
        data = load_untouch_nii(preproc_runs{runnum});
        data.img = data.img(:,:,:,logical(tmask(runIDs==runnum)));
        if runnum == 1
            out = data;
        else
            newtimepoints = size(data.img,4);
            out.img(:,:,:,end+1:end+newtimepoints) = data.img;
        end
        clear data;
    end
    save_untouch_nii(out,[funcvol '_unprocessed_tmasked.nii.gz']);
    
    systemresult = goodvoxels_singlesub([funcvol '_unprocessed_tmasked'], fsLRfolder, goodvoxfolder,subject,sequence_name,systemresult);
    
else
    
    systemresult = goodvoxels_singlesub(funcvol, fsLRfolder, goodvoxfolder,subject,sequence_name,systemresult);
    
end
    
submask = [goodvoxfolder '/' sequence_name '_goodvoxels.nii.gz'];



% Sample volumes to surface, downsample, and smooth
for hem = 1:2
    
    midsurf = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
    midsurf_LR32k = [fsLRfolder '/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    whitesurf = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
    pialsurf = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
    nativedefsphere = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
    outsphere = [fsLRfolder '/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
    surfname = [sequence_name '_' HEMS{hem}];
    
    %disp(['Subject ' subject ': mapping ' HEMS{hem} ' hemisphere data to surface']);
    [systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -volume-to-surface-mapping ' funcvol '.nii.gz ' midsurf ' ' surffuncdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
    
    %disp(['Subject ' subject ': Dilating ' HEMS{hem} ' hemisphere surface timecourse']);
    [systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -metric-dilate ' surffuncdir '/' surfname '.func.gii ' midsurf ' 10 ' surffuncdir '/' surfname '_dil10.func.gii']);
    
    %disp(['Subject ' subject ': Deforming ' HEMS{hem} ' hemisphere timecourse to 32k fs_LR']);
    [systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -metric-resample ' surffuncdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
    
    %disp(['Subject ' subject ': Smoothing ' HEMS{hem} ' hemisphere surface timecourse']);
    [systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii ' num2str(smoothnum) ' ' surffuncdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii']);
    
    surfname_final{hem} = [surffuncdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii'];
    
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['/usr/local/caret/bin_linux64/caret_command -file-convert -format-convert XML_BASE64 ' surfname_final{hem}]);
    
    delete([surffuncdir '/' surfname '.func.gii']);
    delete([surffuncdir '/' surfname '_dil10.func.gii']);
    delete([surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii']);
end



% Smooth data in volume within mask
%disp(['Subject ' subject ': Smoothing functional data within volume mask']);
funcvol_ROIsmooth = [funcvol '_wROI255'];
[systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -volume-smoothing ' funcvol '.nii.gz ' num2str(smoothnum) ' ' funcvol_ROIsmooth '.nii.gz -roi ' subcort_mask]);

delete([funcvol '.nii.gz'])
delete([funcvol '_unprocessed_tmasked.nii.gz'])



% Create cifti timeseries
%disp(['Subject ' subject ': Combining surface and volume data to create cifti timeseries']);
[systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftidir '/' sequence_name '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii -volume ' funcvol_ROIsmooth '.nii.gz ' subcort_mask ' -left-metric ' surfname_final{1} ' -roi-left ' medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' medial_mask_R ' -timestep ' num2str(TR) ' -timestart 0']);

if do_smallwall
    [systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftiswdir '/' sequence_name '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii -volume ' funcvol_ROIsmooth '.nii.gz ' subcort_mask ' -left-metric ' surfname_final{1} ' -roi-left ' sw_medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' sw_medial_mask_R ' -timestep ' num2str(TR) ' -timestart 0']);
end

delete([funcvol_ROIsmooth '.nii.gz'])

