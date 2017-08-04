function systemresult = make_fs_masks(freesurferdir,voxdim,native2MNItransform,template,systemresult)

voxdim = num2str(voxdim);
iter_dilate_wb = 4;
cerebralwm = {'2','41'};
iter_erode_wm = 6;
CSF = {'4','5','14','15','24','43','44'};
iter_erode_CSF = 4;
GM = {'3','42'};

if ~exist('systemresult')
    systemresult = cell(0,2);
end

cd(freesurferdir)
outputdir = [freesurferdir '/nusmask'];
mkdir(outputdir)
[systemresult{end+1,1},systemresult{end+1,2}] = system(['chmod g+w ' outputdir]);

filebase = 'aparc+aseg';
cd([freesurferdir '/mri'])
[systemresult{end+1,1},systemresult{end+1,2}] = system(['mri_convert -rl rawavg.mgz ' filebase '.mgz ' filebase '.nii']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' filebase '.nii ' outputdir]);
cd(outputdir)



%brain mask
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase ' -bin ' filebase '_brainmask']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_brainmask -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_brainmask_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase '_brainmask_mask_' voxdim voxdim voxdim ' -kernel 3D -dilM ' filebase '_brainmask_dil1_mask_' voxdim voxdim voxdim]);
for iter = 2:iter_dilate_wb
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase '_brainmask_dil' num2str(iter-1) '_mask_' voxdim voxdim voxdim ' -kernel 3D -dilM ' filebase '_brainmask_dil' num2str(iter) '_mask_' voxdim voxdim voxdim]);
end


%wm masks
addstr = 'fslmaths ';
for i = 1:length(cerebralwm)
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase ' -thr ' cerebralwm{i} ' -uthr ' cerebralwm{i} ' ' filebase '_reg' cerebralwm{i}]);
    addstr = [addstr filebase '_reg' cerebralwm{i} ' -add '];
end
addstr = [addstr(1:end-5) ' ' filebase '_cerebralwm'];
[systemresult{end+1,1},systemresult{end+1,2}] = system(addstr);
delete([filebase '_reg*.nii.gz'])
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase '_cerebralwm -bin ' filebase '_cerebralwm']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_cerebralwm -ref ' template ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_mask_MNI -interp nearestneighbour']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_cerebralwm -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_ero0_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);

for iter = 1:iter_erode_wm
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase '_cerebralwm -kernel 3D -ero ' filebase '_cerebralwm_ero' num2str(iter)]);
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_cerebralwm_ero' num2str(iter) ' -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_ero' num2str(iter) '_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);
end
%delete([filebase '_cerebralwm.nii.gz'])



%csf masks
addstr = 'fslmaths ';
for i = 1:length(CSF)
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase ' -thr ' CSF{i} ' -uthr ' CSF{i} ' ' filebase '_reg' CSF{i}]);
    addstr = [addstr filebase '_reg' CSF{i} ' -add '];
end
addstr = [addstr(1:end-5) ' ' filebase '_CSF'];
[systemresult{end+1,1},systemresult{end+1,2}] = system(addstr);
delete([filebase '_reg*.nii.gz'])
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase '_CSF -bin ' filebase '_CSF']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_CSF -ref ' template ' -init ' native2MNItransform ' -out ' filebase '_CSF_mask_MNI -interp nearestneighbour']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_CSF -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_CSF_ero0_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);

for iter = 1:iter_erode_wm
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase '_CSF -kernel 3D -ero ' filebase '_CSF_ero' num2str(iter)]);
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_CSF_ero' num2str(iter) ' -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_CSF_ero' num2str(iter) '_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);
end
%delete([filebase '_CSF.nii.gz'])
%delete([filebase '.nii'])



%grey ribbon
filebase = 'aseg';
cd([freesurferdir '/mri'])
[systemresult{end+1,1},systemresult{end+1,2}] = system(['mri_convert -rl rawavg.mgz ' filebase '.mgz ' filebase '.nii']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' filebase '.nii ' outputdir]);
cd(outputdir)

addstr = ['fslmaths '];
for i = 1:length(GM)
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase ' -thr ' GM{i} ' -uthr ' GM{i} ' ' filebase '_reg' GM{i}]);
    addstr = [addstr filebase '_reg' GM{i} ' -add '];
end
addstr = [addstr(1:end-5) ' ' filebase '_GM'];
[systemresult{end+1,1},systemresult{end+1,2}] = system(addstr);
delete([filebase '_reg*.nii.gz'])
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' filebase '_GM -bin ' filebase '_GM']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_GM -ref ' template ' -init ' native2MNItransform ' -out ' filebase '_GM_mask_MNI -interp nearestneighbour']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_GM -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_GM_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);
%delete([filebase '_GM.nii.gz'])
delete([filebase '.nii'])