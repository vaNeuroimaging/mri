%% 120
basedir = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/';
funcdir = [basedir '/FCPROCESS_bandpass_interp_nosmooth'];
outputdir = [funcdir '/surf_timecourses'];
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR';
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
smooth = '2.55';
[subjects tmasks] = textread([basedir '/NEW_nokids_TMASKLIST.txt'],'%s%s');


for s = 1:length(subjects)
    
    disp(['processing subject #' num2str(s) ': ' subjects{s} ])
    for hem = 1:2
              
        midsurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.midthickness.native.surf.gii'];
        midsurf_LR32k = [surfdir '/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        whitesurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.pial.native.surf.gii'];
        nativedefsphere = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        outsphere = [surfdir '/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
   
        
        subfunc = [funcdir '/' subjects{s} '/' subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt'];
        submask = [funcdir '/' subjects{s} '/goodvoxels/' subjects{s} '_goodvoxels.nii.gz'];
        system(['niftigz_4dfp -n ' subfunc ' ' subfunc])
        
        surfname = [subjects{s} '_' HEMS{hem} '_time'];
        disp('Map volume to surface')
        system([workbenchdir '/wb_command -volume-to-surface-mapping ' subfunc '.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
        disp('Dilate surface timecourse')
        system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
        disp('Deform timecourse to 32k fs_LR')
        system([workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
        disp('Smooth surface timecourse')
        system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii ' smooth ' ' outputdir '/' surfname '_dil10_32k_fs_LR_smooth' smooth '.func.gii'])
       % system(['rm ' outputdir '/' surfname '.func.gii'])
       % system(['rm ' outputdir '/' surfname '_dil10.func.gii'])
       % system(['rm ' outputdir '/' surfname '_dil10_smooth' smooth '.func.gii'])
    end
end