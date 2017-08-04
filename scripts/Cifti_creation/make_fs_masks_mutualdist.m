function systemresult = make_fs_masks_mutualdist(freesurferdir,voxdim,native2MNItransform,template,systemresult)
% systemresult = make_fs_masks_mutualdist(freesurferdir,voxdim,native2MNItransform,template,[systemresult])

mutualdist = 6;
voxdim = num2str(voxdim);
iter_dilate_wb = 4;
cerebralwm = {'2','41'};
CSF = {'4','5','14','15','24','43','44'};
GM = {'3','42'};

if ~exist('systemresult')
    systemresult = cell(0,2);
end

cd(freesurferdir)
outputdir = [freesurferdir '/nusmask'];
mkdir(outputdir)
[~,~] = system(['chmod g+w ' outputdir]);

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
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_cerebralwm -ref ' template ' -applyxfm -init ' native2MNItransform ' -out ' filebase '_cerebralwm_mask_MNI -interp nearestneighbour']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_cerebralwm -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);

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
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_CSF -ref ' template ' -applyxfm -init ' native2MNItransform ' -out ' filebase '_CSF_mask_MNI -interp nearestneighbour']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_CSF -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_CSF_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);

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
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_GM -ref ' template ' -applyxfm -init ' native2MNItransform ' -out ' filebase '_GM_mask_MNI -interp nearestneighbour']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_GM -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_GM_mask_' voxdim voxdim voxdim ' -interp nearestneighbour']);

%load Grey and get coords
GM = load_untouch_nii([filebase '_GM_mask_' voxdim voxdim voxdim '.nii.gz']);
xform = [GM.hdr.hist.srow_x; GM.hdr.hist.srow_y; GM.hdr.hist.srow_z];
GMvoxels = find(GM.img);
[GMcoords(:,1),GMcoords(:,2),GMcoords(:,3)] = ind2sub(size(GM.img),GMvoxels);
GMcoords = (xform(:,1:3) * GMcoords'  + repmat(xform(:,4),1,length(GMcoords)))';

filebase = 'aparc+aseg';
%load white and get coords
WM = load_untouch_nii([filebase '_cerebralwm_mask_' voxdim voxdim voxdim '.nii.gz']);
xform = [WM.hdr.hist.srow_x; WM.hdr.hist.srow_y; WM.hdr.hist.srow_z];
WMvoxels = find(WM.img);
[WMcoords(:,1),WMcoords(:,2),WMcoords(:,3)] = ind2sub(size(WM.img),WMvoxels);
WMcoords = (xform(:,1:3) * WMcoords'  + repmat(xform(:,4),1,length(WMcoords)))';

%load csf and get coords
CSF = load_untouch_nii([filebase '_CSF_mask_' voxdim voxdim voxdim '.nii.gz']);
xform = [CSF.hdr.hist.srow_x; CSF.hdr.hist.srow_y; CSF.hdr.hist.srow_z];
CSFvoxels = find(CSF.img);
[CSFcoords(:,1),CSFcoords(:,2),CSFcoords(:,3)] = ind2sub(size(CSF.img),CSFvoxels);
CSFcoords = (xform(:,1:3) * CSFcoords'  + repmat(xform(:,4),1,length(CSFcoords)))';

%find white far from grey and csf
WM2GM_dist = pdist2(WMcoords,GMcoords);
WM2CSF_dist = pdist2(WMcoords,CSFcoords);
invalidWM = any((WM2GM_dist < mutualdist),2) | any((WM2CSF_dist < mutualdist),2);
invalidWMvoxels = WMvoxels(invalidWM);
WM.img(invalidWMvoxels) = 0;
save_untouch_nii(WM,[filebase '_cerebralwm_mask_' voxdim voxdim voxdim '_' num2str(mutualdist) 'mm_fromGMCSF.nii.gz']);
%[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_cerebralwm_' num2str(mutualdist) 'mm_fromGMCSF -ref ' template ' -applyxfm -init ' native2MNItransform ' -out ' filebase '_cerebralwm_mask_MNI_' num2str(mutualdist) 'mm_fromGMCSF -interp nearestneighbour']);
%[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_cerebralwm_' num2str(mutualdist) 'mm_fromGMCSF -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_cerebralwm_mask_' voxdim voxdim voxdim '_' num2str(mutualdist) 'mm_fromGMCSF -interp nearestneighbour']);

%find csf far from grey and white
CSF2GM_dist = pdist2(CSFcoords,GMcoords);
invalidCSF = any((CSF2GM_dist < mutualdist),2) | any((WM2CSF_dist' < mutualdist),2);
invalidCSFvoxels = CSFvoxels(invalidCSF);
CSF.img(invalidCSFvoxels) = 0;
save_untouch_nii(CSF,[filebase '_CSF_mask_' voxdim voxdim voxdim '_' num2str(mutualdist) 'mm_fromGMWM.nii.gz']);
%[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_CSF_' num2str(mutualdist) 'mm_fromGMCSF -ref ' template ' -applyxfm -init ' native2MNItransform ' -out ' filebase '_CSF_mask_MNI_' num2str(mutualdist) 'mm_fromGMCSF -interp nearestneighbour']);
%[systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' filebase '_CSF_' num2str(mutualdist) 'mm_fromGMCSF -ref ' template ' -applyisoxfm ' voxdim ' -init ' native2MNItransform ' -out ' filebase '_CSF_mask_' voxdim voxdim voxdim '_' num2str(mutualdist) 'mm_fromGMCSF -interp nearestneighbour']);






%delete([filebase '_GM.nii.gz'])
%delete([filebase '.nii'])