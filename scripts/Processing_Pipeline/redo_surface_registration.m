subject = 'MAV008';




preprocfolder = ['/home/data/subjects/' subject '/preprocessed/']; mkdir(preprocfolder)
freesurferfolder = ['/home/data/subjects/' subject '/freesurfer/']; 
fsLRfolder = ['/home/data/subjects/' subject '/fs_LR/']; 

T1file = 'T1_biascorr';
template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz';


systemresult = [];
systemresult = make_fs_masks(freesurferfolder,3,[preprocfolder '/T1_2MNI.mat'],template,systemresult);
copyfile([freesurferfolder '/nusmask/aparc+aseg_cerebralwm.nii.gz'],preprocfolder)

mkdir(fsLRfolder)

prevsize = size(systemresult,1);
systemresult = PostFreeSurferPipeline_fsavg2fslr_long_singlesub(subject,T1file,freesurferfolder,fsLRfolder,preprocfolder,systemresult);
if any(cell2mat(systemresult(prevsize+1:end,1)))
    error('Problem detected with surface registration')
    systemresult(prevsize+1:end,2)
end