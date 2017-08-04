function post_fc_processing_batch(post_fc_processing_batch_params_file)
%post_fc_processing_batch(post_fc_processing_batch_params_file)
%
% Run through post fc-processing, creating cifti timeseries from
% unsmoothed fc-processed functional data.
%
% You must have FSLv5.0.6 or later in your path. Earlier versions of FSL
% produce a mysterious alignment error.
%
% Requires the FULL PATH to a parameters file (a .m file) which will be
% executed to load needed parameters, including:
%
% a datalist and a tmasklist, as used by FCPROCESS
% the location of the fc-processed data
% an output folder for the cifti timeseries
% the location of a subcortical mask
% the location of subjects' fs_LR-registered cortical surfaces
% the locations of a priori left and right hemisphere atlas medial wall masks
% the locations of a priori left and right hemisphere "smallwall" medial wall masks
% the sigma of the smoothing kernel to apply to the data
%
%
% Requires the Cifti Resources scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ and subfolders)
%
%EMG 06/24/15



%Find params file
[paramspath,paramsname,paramsextension] = fileparts(post_fc_processing_batch_params_file);
origpath = pwd;
if ~isempty(paramspath)
    cd(paramspath)
end

%Load parameters
params = feval(paramsname);
varnames = fieldnames(params);
for i = 1:length(varnames)
    evalc([varnames{i} ' = params.' varnames{i}]);
end
clear varnames params

cd(origpath)

do_smallwall = true;
if (~exist('sw_medial_mask_L','var')) || (~exist('sw_medial_mask_R','var')) || isempty(sw_medial_mask_L) || isempty(sw_medial_mask_R)
    do_smallwall = false;
end
    


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
[~, tmasks] = textread(tmasklist,'%s %s');

prevstring = [];

delete(outputdatalistname);
fid = fopen(outputdatalistname,'at'); %open the output file for writing
fclose(fid);

if do_smallwall
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
        if do_smallwall
            mkdir(ciftiswdir);
        end
    end
    
    
    %Convert and copy data over
    string = ['Subject ' subject ': converting and copying fc-processed functional data'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    subfunc = [fcprocessed_funcdata_dir '/' subject '/' subject fcprocessed_funcdata_dir_suffix];
    if ~exist([subfunc '.nii.gz'])
        error(['fc-processed subject data ' subfunc ' does not exist!'])
    end
    funcvol = [surffuncdir '/' subject '_funcvol'];
    copyfile(subfunc,funcvol)    
    
    %Remove NaNs from the data
    evalc(['!fslmaths ' funcvol ' -nan ' funcvol]);
        
    
    %Identify high-SNR voxels
    string = ['Subject ' subject ': calculating high-SNR voxels from unprocessed functional data'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    [ign, runs] = textread(prmfiles{s},'%s%s','delimiter','(');
    runs = str2num(runs{1}(1:end-1));
    fslstr = ['!fslmerge -t ' funcvol '_unprocessed'];
    for runnum = 1:length(runs)
        %unproc_rundata = [funcdirs{s} '/' subject '/bold' num2str(runnum) '/' subject '_b' num2str(runnum) '_xr3d_333.4dfp.img'];
        unproc_rundata = ['/data/subjects/' subject '/preprocessed/' runs{runnum}];
        if ~exist(unproc_rundata)
            error(['Unprocessed subject data ' unproc_rundata ' does not exist!'])
        end
        fslstr = [fslstr ' ' unproc_rundata];
    end
    evalc(fslstr);
       
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
    
    if do_smallwall
        evalc(['!' workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftiswdir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii -volume ' funcvol_ROIsmooth '.nii.gz ' subcort_mask ' -left-metric ' surfname_final{1} ' -roi-left ' sw_medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' sw_medial_mask_R ' -timestep ' num2str(TR) ' -timestart 0']);
        dlmwrite(outputdatalistname_sw,[subject '  ' ciftiswdir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii'],'-append','delimiter','');
    end
    
    delete([funcvol_ROIsmooth '.nii.gz'])
    
    
    
    
end
disp(' ')