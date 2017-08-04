function post_fc_processing_batch_v2(post_fc_processing_batch_params_file)
%post_fc_processing_batch(post_fc_processing_batch_params_file)
%
% Run through post fc-processing, creating cifti timeseries from
% unsmoothed fc-processed functional data.
%
% Requires a parameters file (a .m file) which will be executed to load
% needed parameters, including:
%
% a datalist and a tmasklist, as used by FCPROCESS
% the location of the fc-processed data
% an output folder for the cifti timeseries
% the location of a subcortical mask
% the location of subjects' fs_LR-registered cortical surfaces
% the locations of a priori left and right hemisphere medial wall masks
% the sigma of the smoothing kernel to apply to the data
%
%EMG 06/24/15



%Load parameters
[paramspath,paramsname,paramsextension] = fileparts(post_fc_processing_batch_params_file);
origpath = pwd;
cd(paramspath)

[datalist,tmasklist,fcprocessed_funcdata_dir,fcprocessed_funcdata_dir_suffix,outfolder,subcort_mask,fs_LR_surfdir,medial_mask_L,medial_mask_R,sw_medial_mask_L,sw_medial_mask_R,smoothnum] = feval(paramsname);

cd(origpath)



%-------------------------------------------------------------------------

warning off
mkdir(outfolder)
surffuncdir = [outfolder '/surf_timecourses/'];
ciftidir = [outfolder '/cifti_timeseries_normalwall/'];
ciftiswdir = [outfolder '/cifti_timeseries_smallwall/'];
goodvoxfolder = [outfolder '/goodvoxels/'];
outputdatalistname = [outfolder '/cifti_datalist.txt'];
outputdatalistname_sw = [outfolder '/cifti_datalist_smallwall.txt'];

workbenchdir = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/workbench/bin_linux64/';
HEMS = {'L';'R'};


[funcdirs, subjects, prmfiles, TRs, skip] = textread(datalist,'%s %s %s %s %s');
[ign, tmasks] = textread(tmasklist,'%s %s');

prevstring = [];

delete(outputdatalistname);
fid = fopen(outputdatalistname,'at'); %open the output file for writing
fclose(fid);

if ~isempty(sw_medial_mask_L)
    delete(outputdatalistname_sw);
    fid = fopen(outputdatalistname_sw,'at'); %open the output file for writing
    fclose(fid);
end


for s = 1:length(subjects)
    subject = subjects{s};
    tmask = load(tmasks{s});
    TR = TRs{s};
    
    
    %Make folders
    if s==1;
        mkdir(surffuncdir);
        mkdir(ciftidir);
        mkdir(goodvoxfolder);
        if ~isempty(sw_medial_mask_L) 
            mkdir(ciftiswdir);
        end
    end
    
    
    %Convert and copy data over
    string = ['Subject ' subject ': converting and copying fc-processed functional data'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    subfunc = [fcprocessed_funcdata_dir '/' subject '/' subject fcprocessed_funcdata_dir_suffix];
    if ~exist([subfunc '.4dfp.img'])
        error(['fc-processed subject data ' subfunc ' does not exist!'])
    end
    funcvol = [surffuncdir '/' subject '_funcvol'];
    evalc(['!niftigz_4dfp -n ' subfunc ' ' funcvol]);
    
    
    
    %Identify high-SNR voxels
    string = ['Subject ' subject ': calculating high-SNR voxels from unprocessed functional data'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    [ign, runs] = textread(prmfiles{s},'%s%s','delimiter','(');
    runs = str2num(runs{1}(1:end-1));
    fslstr = ['!fslmerge -t ' funcvol '_unprocessed'];
    for runnum = runs
        unproc_rundata = [funcdirs{s} '/' subject '/bold' num2str(runnum) '/' subject '_b' num2str(runnum) '_xr3d_333.4dfp.img'];
        if ~exist(unproc_rundata)
            error(['Unprocessed subject data ' unproc_rundata ' does not exist!'])
        end
        evalc(['!niftigz_4dfp -n ' unproc_rundata ' ' surffuncdir '/temp' num2str(runnum)]);
        fslstr = [fslstr ' ' surffuncdir '/temp' num2str(runnum) '.nii.gz'];
    end
    evalc(fslstr);
    delete([surffuncdir '/temp*.nii.gz'])
       
    goodvoxels([funcvol '_unprocessed'], tmask, fs_LR_surfdir, goodvoxfolder,subject);
    
    submask = [goodvoxfolder '/' subject '_goodvoxels.nii.gz'];
        
    
    
    % Sample volumes to surface, downsample, and smooth
    for hem = 1:2
        
        midsurf = [fs_LR_surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
        midsurf_LR32k = [fs_LR_surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        whitesurf = [fs_LR_surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [fs_LR_surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
        nativedefsphere = [fs_LR_surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        outsphere = [fs_LR_surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        surfname = [subject '_' HEMS{hem}];
        
        string = ['Subject ' subject ': mapping ' HEMS{hem} ' hemisphere data to surface'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        evalc(['!' workbenchdir '/wb_command -volume-to-surface-mapping ' funcvol '.nii.gz ' midsurf ' ' surffuncdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
                
        string = ['Subject ' subject ': Dilating ' HEMS{hem} ' hemisphere surface timecourse'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        evalc(['!' workbenchdir '/wb_command -metric-dilate ' surffuncdir '/' surfname '.func.gii ' midsurf ' 10 ' surffuncdir '/' surfname '_dil10.func.gii']);
        
        string = ['Subject ' subject ': Deforming ' HEMS{hem} ' hemisphere timecourse to 32k fs_LR'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        evalc(['!' workbenchdir '/wb_command -metric-resample ' surffuncdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
        
        string = ['Subject ' subject ': Smoothing ' HEMS{hem} ' hemisphere surface timecourse'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        evalc(['!' workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii ' num2str(smoothnum) ' ' surffuncdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii']);
        
        surfname_final{hem} = [surffuncdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii'];
                
        evalc(['!caret_command64 -file-convert -format-convert XML_BASE64 ' surfname_final{hem}]);
        
        delete([surffuncdir '/' surfname '.func.gii']);
        delete([surffuncdir '/' surfname '_dil10.func.gii']);
        delete([surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii']);
    end
    
    
    
    % Smooth data in volume within mask
    string = ['Subject ' subject ': Smoothing functional data within volume mask'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    funcvol_ROIsmooth = [funcvol '_wROI255'];
    evalc(['!' workbenchdir '/wb_command -volume-smoothing ' funcvol '.nii.gz ' num2str(smoothnum) ' ' funcvol_ROIsmooth '.nii.gz -roi ' subcort_mask]);
    
    delete([funcvol '.nii.gz'])
    delete([funcvol '_unprocessed.nii.gz'])
    
    
    
    % Create cifti timeseries
    string = ['Subject ' subject ': Combining surface and volume data to create cifti timeseries'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    evalc(['!' workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftidir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii -volume ' funcvol_ROIsmooth '.nii.gz ' subcort_mask ' -left-metric ' surfname_final{1} ' -roi-left ' medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' medial_mask_R ' -timestep ' num2str(TR) ' -timestart 0']);
    dlmwrite(outputdatalistname,[subject '  ' ciftidir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii'],'-append','delimiter','');
    
    if ~isempty(sw_medial_mask_L) 
        evalc(['!' workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftiswdir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii -volume ' funcvol_ROIsmooth '.nii.gz ' subcort_mask ' -left-metric ' surfname_final{1} ' -roi-left ' sw_medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' sw_medial_mask_R ' -timestep ' num2str(TR) ' -timestart 0']);
        dlmwrite(outputdatalistname_sw,[subject '  ' ciftiswdir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii'],'-append','delimiter','');
    end
    
    delete([funcvol_ROIsmooth '.nii.gz'])
    
    
    
    
end
disp(' ')