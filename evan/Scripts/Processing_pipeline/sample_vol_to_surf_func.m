function sample_vol_to_surf_func(subject,funcdir,goodvoxfolder,smoothnum,outputdir,surfdir)
% Function for sampling volumes to the surface following 4 steps:
% 1. Volume-to-surface mapping using ribbon-constrained sampling to native
% surface
% 2. Dilation of data on surface
% 3. Downsample to 32k surface
% 4. Smooth data along surface
%
% If a subject is specified as the last variable then tmasklist is treated
% as different sessions of a single subject, rather than different subjects
% TOL, 09/14

workbenchdir = '/data/cn4/evan/workbench/bin_linux64/';
HEMS = {'L';'R'};

system(['mkdir ' outputdir])



session = ['bold' num2str(boldruns(s))];
disp(['processing subject #' subject])
subfunc = [funcdir '/' subject '/' subject '_333_zmdt_resid_bpss_zmdt_g7'];
if ~isempty(goodvoxfolder)
    submask = [goodvoxfolder '/' session '_goodvoxels.nii.gz'];
end
%submask = [bolddir '/' session '/goodvoxels/' session '_goodvoxels.nii.gz'];
system(['niftigz_4dfp -n ' subfunc ' ' outputdir '/' subject '_funcvol']);

for hem = 1:2
    
    midsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
    midsurf_LR32k = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    whitesurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
    pialsurf = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
    nativedefsphere = [surfdir '/' subject '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
    outsphere = [surfdir '/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
    
    surfname = [subject '_' HEMS{hem}];
    disp('Map volume to surface')
    if ~isempty(goodvoxfolder)
        system([ workbenchdir '/wb_command -volume-to-surface-mapping ' outputdir '/' subject '_funcvol.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
    else
        system([ workbenchdir '/wb_command -volume-to-surface-mapping ' outputdir '/' subject '_funcvol.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf]);
    end
    disp('Dilate surface timecourse')
    system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
    disp('Deform timecourse to 32k fs_LR')
    system([ workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
    disp('Smooth surface timecourse')
    system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii ' num2str(smoothnum) ' ' outputdir '/' surfname '_dil10_32k_fs_LR_smooth' num2str(smoothnum) '.func.gii'])
    
    system(['rm ' outputdir '/' surfname '.func.gii']);
    system(['rm ' outputdir '/' surfname '_dil10.func.gii']);
    system(['rm ' outputdir '/' surfname '_dil10_32k_fs_LR.func.gii']);
end



    
