% Run through fc post-processing, preparing for parcellation and infomap
TR = 2.2;
smoothnum = 2.55;
basedir = '/net/nil-bluearc/GMT/Laumann/MSC/';
%basedir = '/data/heisenberg/data1/MSC';
subjects = {'MSC01';'MSC02';'MSC03';'MSC04';'MSC05';'MSC06'};

for s  = 1:length(subjects)
    subject = subjects{s};
    bolddir = [basedir '/' subject '/Functionals'];
    funcdir = [bolddir '/FCPROCESS_SCRUBBED_UWRPMEAN'];
    maskdir = [basedir '/' subject '/subcortical_mask_atlas_freesurf'];
    freesurfdir = ['/net/nil-bluearc/GMT/Laumann/MSC/freesurfer_native/' subject];
    surfdir = '/net/nil-bluearc/GMT/Laumann/MSC/freesurfer_native/FREESURFER_fs_LR';
    %freesurfdir = ['/data/heisenberg/data1/MSC/fs5.3/' subject 'edits'];
    %surfdir = ['/data/heisenberg/data1/MSC/fs5.3/FREESURFER_fs_LR'];
    tmasklist = [basedir '/' subject '/' subject '_TMASKLIST.txt'];
    
    % Create subcortical mask
    %      mkdir(maskdir)
    %      t4file = [basedir '/' subject '/T1/' subject '_mpr1T_debias_to_TRIO_Y_NDC_noscale_t4'];
    %      create_subcortical_mask_func_MSC(freesurfdir,t4file,maskdir)
    
    create_subcortical_mask_func(freesurfdir,maskdir)
    
    % Sample volumes to surface, downsample, and smooth
    surffuncdir = [funcdir '/surf_timecourses_atlas_freesurf'];
    sample_vol_to_surf_func_MSC(tmasklist,funcdir,bolddir,smoothnum,surffuncdir,surfdir,subject)
    
    % Smooth data in volume within mask
    fcprocess_smooth_volume_wROI_func(tmasklist,funcdir,maskdir,smoothnum)
    
    % Create medial wall mask based on surface projection
    create_medialmask_from_proj_func(tmasklist,funcdir,maskdir)
    
    % Create cifti timeseries, normalwall
    outputdir = [funcdir '/cifti_timeseries_normalwall_atlas_freesurf'];
    create_cifti_timeseries_manysubs_func(tmasklist,funcdir,surffuncdir,outputdir,TR,'normalwall',maskdir)
    
    % Create cifti timeseries, smallwall
    outputdir = [funcdir '/cifti_timeseries_smallwall_atlas_freesurf'];
    create_cifti_timeseries_manysubs_func(tmasklist,funcdir,surffuncdir,outputdir,TR,'smallwall',maskdir)
end