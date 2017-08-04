subjectlist = '/data/cn4/evan/Task_parcellation/HCPSubjects/SubjectList3.txt';

subjects = textread(subjectlist,'%s');

betalist = '/data/cn4/evan/Task_parcellation/HCPSubjects/MainTaskContrasts_torun.txt';

[ign betanums ign ign ign ign] = textread(betalist,'%s%s%s%s%s%s');

for i=2:length(betanums)
    betastouse(i-1) = str2num(betanums{i});
end

datafolder = '/data/cn4/evan/Task_parcellation/HCPSubjects/';

basesurfdir = '/data/cn4/evan/HCP_Structural/';

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';


thisdir = pwd;

cd(datafolder)

for s = 1:length(subjects)
    
%     surfdir = [basesurfdir subjects{s} '/MNINonLinear/'];
%     
     gunzip([subjects{s} '_mergeFile.nii.gz']);
%     
%     midsurf = [surfdir '/Native/' subjects{s} '.L.midthickness.native.surf.gii'];
%     whitesurf = [surfdir '/Native/' subjects{s} '.L.white.native.surf.gii'];
%     pialsurf = [surfdir '/Native/' subjects{s} '.L.pial.native.surf.gii'];
%     
%     system([workbenchdir '/wb_command -volume-to-surface-mapping ' subjects{s} '_mergeFile.nii ' midsurf ' ' subjects{s} '_mergeFile_L.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -voxel-subdiv 5']);
%     
%     % Deform from native to 32K fs_LR surface
%     cd([surfdir '/fsaverage_LR32k/'])
%     
%     system(['caret_command64 -deformation-map-apply native232k_fs_LR.L.deform_map METRIC_AVERAGE_TILE ' datafolder '/' subjects{s} '_mergeFile_L.func.gii ' datafolder '/' subjects{s} '_mergeFile_L_32k_fsLR.func.gii']);
    
    cd(datafolder)

%     system(['caret_command64 -file-convert -format-convert ASCII ' subjects{s} '_mergeFile_L_32k_fsLR.func.gii'])
%     system(['awk ''NF > 25'' ' subjects{s} '_mergeFile_L_32k_fsLR.func.gii > ' subjects{s} '_mergeFile_L_32k_fsLR_noHEAD.func.gii'])
%     tempdata = load([subjects{s} '_mergeFile_L_32k_fsLR_noHEAD.func.gii']);
%     data(:,:,s) = tempdata(:,betastouse+1);
%     
%     system(['caret_command64 -file-convert -format-convert XML ' subjects{s} '_mergeFile_L_32k_fsLR.func.gii'])
    
    tempdata = load_nifti([subjects{s} '_mergeFile.nii']);
    
    voldata(:,:,:,:,s) = tempdata.vol(:,:,:,betastouse);


    delete([subjects{s} '_mergeFile.nii']);
    delete([subjects{s} '_mergeFile_L.func.gii']);
    delete([subjects{s} '_mergeFile_L_32k_fsLR_noHEAD.func.gii']);
    
end

%  Meanbetas = squeeze(mean(data,3));
%  
%  save(gifti(Meanbetas),'Meansub_allbetas_L.func.gii')
%  
 Meanbetasvol = squeeze(mean(voldata,5));
 tempdata.vol = Meanbetasvol;
 save_nifti(tempdata,'Meansub_allbetas.nii');
 

cd(thisdir)