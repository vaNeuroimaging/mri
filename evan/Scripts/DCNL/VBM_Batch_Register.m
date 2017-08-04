

%
%Created by E. Gordon 5/08/08


warning off

subjects = {'CON_016.nii' 'CON_017.nii' 'CON_018.nii' 'CON_019.nii' 'CON_020.nii' 'CON_021.nii' 'CON_022.nii' 'CON_023.nii' 'CON_024.nii' 'CON_025.nii' 'CON_026.nii' 'CON_027.nii' 'CON_028.nii' 'CON_029.nii' 'CON_030.nii'};
%subjects = {'ASD_001.nii' 'ASD_002.nii' 'ASD_003.nii' 'ASD_004.nii' 'ASD_005.nii' 'ASD_006.nii' 'ASD_007.nii' 'ASD_008.nii' 'ASD_009.nii' 'ASD_010.nii' 'ASD_011.nii' 'ASD_012.nii' 'ASD_013.nii' 'ASD_014.nii' 'ASD_015.nii' 'ASD_016.nii' 'ASD_017.nii' 'ASD_018.nii' 'ASD_019.nii' 'ASD_020.nii' 'ASD_021.nii' 'ASD_022.nii' 'ASD_023.nii' 'ASD_024.nii' 'ASD_025.nii' 'ASD_026.nii' 'ASD_027.nii' 'ASD_028.nii' 'ASD_029.nii' 'ASD_030.nii' 'CON_001.nii' 'CON_002.nii' 'CON_003.nii' 'CON_004.nii' 'CON_005.nii' 'CON_006.nii' 'CON_007.nii' 'CON_008.nii' 'CON_009.nii' 'CON_010.nii' 'CON_011.nii' 'CON_012.nii' 'CON_013.nii' 'CON_014.nii' 'CON_015.nii' 'CON_016.nii' 'CON_017.nii' 'CON_018.nii' 'CON_019.nii' 'CON_020.nii' 'CON_021.nii' 'CON_022.nii' 'CON_023.nii' 'CON_024.nii' 'CON_025.nii' 'CON_026.nii' 'CON_027.nii' 'CON_028.nii' 'CON_029.nii' 'CON_030.nii'};

for subject = 1:length(subjects)
    load BATCH_register.mat;
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = ['/fmri/data3/Devon/VBM/ADHD+ASD/DARTEL/SPM8_rerun/' subjects{subject} ',1'];
    save register matlabbatch;
    spm_jobman('run','register.mat');
    
end

