function make_fs_masks(freesurferdir,voxdim,native2MNItransform,template,systemresult)

voxdim = num2str(voxdim);

cerebralwm = {'2','41'};
iter_erode_wm = 6;
CSF = {'4','5','14','15','24','43','44'};
iter_erode_CSF = 4;
GM = {'3','42'};

cd(freesurferdir)
outputdir = [freesurferdir '/nusmask'];
mkdir(outputdir)
system(['chmod g+w ' outputdir])

filebase = 'aparc+aseg';
cd([freesurferdir '/mri'])
system(['mri_convert -rl rawavg.mgz ' filebase '.mgz ' filebase '.nii'])
system(['mv ' filebase '.nii ' outputdir])
cd(outputdir)



%brain mask
system(['fslmaths ' filebase ' -bin ' filebase '_brainmask'])
system(['flirt -in ' filebase '_brainmask -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_brainmask_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);



%wm masks
addstr = ['fslmaths '];
for i = 1:length(cerebralwm)
    system(['fslmaths ' filebase ' -thr ' cerebralwm{i} ' -uthr ' cerebralwm{i} ' ' filebase '_reg' cerebralwm{i}])
    addstr = [addstr filebase '_reg' cerebralwm{i} ' -add '];
end
addstr = [addstr(1:end-5) ' ' filebase '_cerebralwm'];
system(addstr)
delete([filebase '_reg*.nii.gz'])
system(['fslmaths ' filebase '_cerebralwm -bin ' filebase '_cerebralwm'])
system(['flirt -in ' filebase '_cerebralwm -ref ' template ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_mask_MNI -interp nearestneighbour']);
system(['flirt -in ' filebase '_cerebralwm -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_ero0_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);

for iter = 1:iter_erode_wm
    system(['fslmaths ' filebase '_cerebralwm -kernel 3D -ero ' filebase '_cerebralwm_ero' num2str(iter)])
    system(['flirt -in ' filebase '_cerebralwm_ero' num2str(iter) ' -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_ero' num2str(iter) '_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);
end
%delete([filebase '_cerebralwm.nii.gz'])



%csf masks
addstr = ['fslmaths '];
for i = 1:length(CSF)
    system(['fslmaths ' filebase ' -thr ' CSF{i} ' -uthr ' CSF{i} ' ' filebase '_reg' CSF{i}])
    addstr = [addstr filebase '_reg' CSF{i} ' -add '];
end
addstr = [addstr(1:end-5) ' ' filebase '_CSF'];
system(addstr)
delete([filebase '_reg*.nii.gz'])
system(['fslmaths ' filebase '_CSF -bin ' filebase '_CSF'])
system(['flirt -in ' filebase '_CSF -ref ' template ' -init ' native2MNItransform ' -out ' filebase '_CSF_mask_MNI -interp nearestneighbour']);
system(['flirt -in ' filebase '_CSF -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_CSF_ero0_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);

for iter = 1:iter_erode_wm
    system(['fslmaths ' filebase '_CSF -kernel 3D -ero ' filebase '_CSF_ero' num2str(iter)])
    system(['flirt -in ' filebase '_CSF_ero' num2str(iter) ' -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_CSF_ero' num2str(iter) '_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);
end
%delete([filebase '_CSF.nii.gz'])
delete([filebase '.nii'])



%grey ribbon
filebase = 'aseg';
cd([freesurferdir '/mri'])
system(['mri_convert -rl rawavg.mgz ' filebase '.mgz ' filebase '.nii'])
system(['mv ' filebase '.nii ' outputdir])
cd(outputdir)

addstr = ['fslmaths '];
for i = 1:length(GM)
    system(['fslmaths ' filebase ' -thr ' GM{i} ' -uthr ' GM{i} ' ' filebase '_reg' GM{i}])
    addstr = [addstr filebase '_reg' GM{i} ' -add '];
end
addstr = [addstr(1:end-5) ' ' filebase '_GM'];
system(addstr)
delete([filebase '_reg*.nii.gz'])
system(['fslmaths ' filebase '_GM -bin ' filebase '_GM'])
system(['flirt -in ' filebase '_GM -ref ' template ' -init ' native2MNItransform ' -out ' filebase '_GM_mask_MNI -interp nearestneighbour']);
system(['flirt -in ' filebase '_GM -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_GM_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);
%delete([filebase '_GM.nii.gz'])
delete([filebase '.nii'])