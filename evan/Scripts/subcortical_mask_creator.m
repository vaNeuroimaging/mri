outfolder = '/data/cn5/selfRegulation/V4Process_nosmooth/';

[subjects ign] = textread(['/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120.txt'],'%s%s');

subcortical = read_4dfpimg('/data/cn4/evan/ROIs/mode_subcortical_label_LR_333.4dfp.img');
subcortical_labels = unique(subcortical); subcortical_labels(subcortical_labels==0) = [];

for s = 1:length(subjects)
    subname = subjects{s};
    disp(subname)
    segfilefolder = ['/data/cn4/segmentation/freesurfer5_supercomputer/' subname '/mri/'];
    
    system(['mri_convert -rl ' segfilefolder 'rawavg.mgz ' segfilefolder 'aparc.a2009s+aseg.mgz ' outfolder '/Temp.nii'])
    system(['nifti_4dfp -4 ' outfolder '/Temp.nii ' outfolder '/Temp.4dfp.img'])
    system(['t4img_4dfp none ' outfolder 'Temp.4dfp.img ' outfolder 'Temp_333.4dfp.img -0333'])
    
    
    wmparc(:,s) = read_4dfpimg([outfolder 'Temp_333.4dfp.img']);
    [voxelsize frames etype] = read_4dfpifh([outfolder 'Temp_333.4dfp.ifh']);
    
    
end

%%
mode_wmparc = mode(wmparc,2);

%%

mode_subcortical = zeros(size(mode_wmparc));
for labelnum = 1:length(subcortical_labels)
    mode_subcortical(mode_wmparc==subcortical_labels(labelnum)) = subcortical_labels(labelnum);
end

write_4dfpimg(mode_subcortical,[outfolder 'mode_subcortical_LR_333.4dfp.img'],etype)
write_4dfpifh([outfolder 'mode_subcortical_LR_333.4dfp.ifh'],1,etype)

system(['nifti_4dfp -n ' outfolder 'mode_subcortical_LR_333.4dfp.img ' outfolder 'mode_subcortical_LR_333.nii'])
system(['wb_command -volume-label-import ' outfolder 'mode_subcortical_LR_333.nii /data/cn4/laumannt/subcortical_mask/FreeSurferSubcorticalLabelTableLut_nobrainstem_LR.txt ' outfolder 'mode_subcortical_label_LR_333.nii'])
