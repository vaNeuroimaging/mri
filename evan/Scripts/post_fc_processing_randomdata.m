% Run through fc post-processing, preparing for parcellation and infomap
warning off

%SPECIFY PARAMETERS

smoothnum = 2.55;

datalist = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';

outfolder = '/data/cn4/evan/Temp/surface_parcellation_randomdata_innative/';

surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/';

medial_mask_L = '/data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii';
medial_mask_R = '/data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii';

ntimepoints = 200;

%-------------------------------------------------------------------------


surffuncdir = [outfolder '/surf_timecourses/'];
mkdir(surffuncdir)
ciftidir = [outfolder '/cifti_timeseries/'];
mkdir(ciftidir)

workbenchdir = '/data/cn4/evan/workbench/bin_linux64/';
HEMS = {'L';'R'};

[subjects ign] = textread(datalist,'%s %s');
subjects = subjects(1:120);
%[funcdirs, subjects, prmfiles, TRs, skip] = textread(datalist,'%s %s %s %s %s');
%[ign, tmasks] = textread(tmasklist,'%s %s');

prevstring = [];
TR = 2.5;

for s  = 1:length(subjects)
    subject = subjects{s};
    
    
    % Sample volumes to surface, downsample, and smooth
    for hem = 1%1:2
        
        midsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
        midsurf_LR32k = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        whitesurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
        nativedefsphere = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        outsphere = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        surfname = [subject '_' HEMS{hem}];
        
        string = ['Subject ' subject ': Making random ' HEMS{hem} ' hemisphere timecourse'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        surf = gifti(midsurf);
        data = rand(length(surf.vertices),ntimepoints);
        save(gifti(single(data)),[surffuncdir '/' surfname '.func.gii'],'ExternalFileBinary');
        
%         string = ['Subject ' subject ': mapping ' HEMS{hem} ' hemisphere data to surface'];
%         fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%         prevstring = string;
%         evalc(['!' workbenchdir '/wb_command -volume-to-surface-mapping ' funcvol '.nii.gz ' midsurf ' ' surffuncdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
%                 
%         string = ['Subject ' subject ': Dilating ' HEMS{hem} ' hemisphere surface timecourse'];
%         fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%         prevstring = string;
%         evalc(['!' workbenchdir '/wb_command -metric-dilate ' surffuncdir '/' surfname '.func.gii ' midsurf ' 10 ' surffuncdir '/' surfname '_dil10.func.gii']);
%         
%         string = ['Subject ' subject ': Deforming ' HEMS{hem} ' hemisphere timecourse to 32k fs_LR'];
%         fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%         prevstring = string;
%         evalc(['!' workbenchdir '/wb_command -metric-resample ' surffuncdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
%         
        string = ['Subject ' subject ': Smoothing ' HEMS{hem} ' hemisphere surface timecourse'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        evalc(['!' workbenchdir '/wb_command -metric-smoothing ' midsurf ' ' surffuncdir '/' surfname '.func.gii ' num2str(smoothnum) ' ' surffuncdir '/' surfname '_smooth' num2str(smoothnum) '.func.gii']);
        
        %surfname_final{hem} = [surffuncdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii'];
        
        surfname_final{hem} = [surffuncdir '/' surfname '_smooth' num2str(smoothnum) '.func.gii'];
                
        evalc(['!caret_command64 -file-convert -format-convert XML_BASE64 ' surfname_final{hem}]);
        
        delete([surffuncdir '/' surfname '.func.gii']);
        %delete([surffuncdir '/' surfname '_dil10.func.gii']);
        delete([surffuncdir '/' surfname '.func.dat']);
        delete([surffuncdir '/' surfname '_dil10_32k_fs_LR.func.gii']);
    end
    
    
    
%     % Smooth data in volume within mask
%     string = ['Subject ' subject ': Smoothing functional data within volume mask'];
%     fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%     prevstring = string;
%     funcvol_ROIsmooth = [funcvol '_wROI255'];
%     evalc(['!' workbenchdir '/wb_command -volume-smoothing ' funcvol '.nii.gz ' num2str(smoothnum) ' ' funcvol_ROIsmooth '.nii.gz -roi ' subcort_maskdir '/subcortical_mask_LR_333.nii']);
%     delete([funcvol '.nii.gz'])
%     delete([funcvol '_unprocessed.nii.gz'])
    
    
    
    % Create cifti timeseries
    string = ['Subject ' subject ': Combining surface data to create cifti timeseries'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    %evalc(['!' workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftidir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii -volume ' funcvol_ROIsmooth '.nii.gz ' subcort_maskdir '/subcortical_mask_LR_333.nii -left-metric ' surfname_final{1} ' -roi-left ' medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' medial_mask_R ' -timestep ' num2str(TR) ' -timestart 0']);
    %evalc(['!' workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftidir '/' subject '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(smoothnum) '.dtseries.nii -left-metric ' surfname_final{1} ' -roi-left ' medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' medial_mask_R ' -timestep ' num2str(TR) ' -timestart 0']);
    
    %evalc(['!' workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftidir '/' subject '_LR_surf_smooth' num2str(smoothnum) '.dtseries.nii -left-metric ' surfname_final{1} ' -right-metric ' surfname_final{2} ' -timestep ' num2str(TR) ' -timestart 0']);
    evalc(['!' workbenchdir '/wb_command -cifti-create-dense-timeseries ' ciftidir '/' subject '_L_surf_smooth' num2str(smoothnum) '.dtseries.nii -left-metric ' surfname_final{1} ' -timestep ' num2str(TR) ' -timestart 0']);
    %delete([funcvol_ROIsmooth '.nii.gz'])
    delete([surffuncdir '/' subject '_L_smooth' num2str(smoothnum) '.func.gii'])
    delete([surffuncdir '/' subject '_R_smooth' num2str(smoothnum) '.func.gii'])
    
end