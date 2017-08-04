%%
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
        whitesurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [surfdir '/' subjects{s} '/7112b_fs_LR/Native/' subjects{s} '.' HEMS{hem} '.pial.native.surf.gii'];
        
        subfunc = [funcdir '/' subjects{s} '/' subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt'];
        submask = [funcdir '/' subjects{s} '/goodvoxels/' subjects{s} '_goodvoxels.nii.gz'];
        system(['niftigz_4dfp -n ' subfunc ' ' subfunc])
        
        surfname = [subjects{s} '_' HEMS{hem} '_time'];
        disp('Map volume to surface')
        system([workbenchdir '/wb_command -volume-to-surface-mapping ' subfunc '.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
        disp('Dilate surface timecourse')
        system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
        disp('Smooth surface timecourse')
        system([workbenchdir '/wb_command -metric-smoothing ' midsurf ' ' outputdir '/' surfname '_dil10.func.gii ' smooth ' ' outputdir '/' surfname '_dil10_smooth' smooth '.func.gii'])
        
        disp('Deform timecourse to 32k fs_LR')
        cd([surfdir '/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/'])
        system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' outputdir '/' surfname '_dil10_smooth' smooth '.func.gii ' outputdir '/' surfname '_dil10_smooth' smooth '_32k_fsLR.func.gii']);
        system(['rm ' outputdir '/' surfname '.func.gii'])
        system(['rm ' outputdir '/' surfname '_dil10.func.gii'])
        system(['rm ' outputdir '/' surfname '_dil10_smooth' smooth '.func.gii'])
    end
end
