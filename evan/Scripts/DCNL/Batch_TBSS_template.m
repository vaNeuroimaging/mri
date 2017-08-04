%Batch_TBSS
%
%Conducts tract-based spatial statistics on DTI data for a given group of
%subjects.  Splits resulting file into individual subject files and
%converts to Nifti for easy loading into SPM group analysis
%
%Subjects and paths are specified at the top of the script.  IMPORTANT:
%Subjects must be listed in alphabetical order!
%
%Requires FSL to be loaded on your computer.
%
%
%Created by E. Gordon 02/11



%USER INPUT (and more below)
%--------------------------------------------------------------------------

%Subjects (MUST be in alphabetical order!!!)
subjects = {'101','102','112','113','118','120','122','125','126','127','132','133','137','138','147','150','151','154','156','159','160','161','162','172','181','182','187','202','207','208','211','214','215','221','222','225','227','229','232','233','242','250','251','253','254','255','256','258','260','261','264','270','272','274','279','281','283','292','301','307','322','327'};

outputfolder = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/';

%END USER INPUT
%--------------------------------------------------------------------------

try rmdir(outputfolder,'s'); catch; end
mkdir(outputfolder);


for subject = 1:length(subjects)
    
    subjid = subjects{subject};
    
    %USER INPUT
    %--------------------------------------------------------------------------
    
    %Location of this subject's FA image
    subjectFA = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/bothruns/dti_FA.nii.gz'];
    
    %END USER INPUT
    %--------------------------------------------------------------------------    
    
    copyfile(subjectFA, [outputfolder '/' subjid '_dti_FA.nii.gz']);
    
end

cd(outputfolder);

disp('Step 1: Preprocessing')

eval('!tbss_1_preproc *.nii.gz');

disp('Step 2: Registration')
 
eval('!tbss_2_reg -T');

disp('Step 3: Skeletonize mean image')
 
eval('!tbss_3_postreg -S');

disp('Step 4: Skeletonize individual images')
 
eval('!tbss_4_prestats .2');


copyfile([outputfolder 'stats/all_FA_skeletonised.nii.gz'],[outputfolder 'stats/indsubs_FA_skeletonized.nii.gz']);

eval(['!fslsplit ' outputfolder '/stats/indsubs_FA_skeletonized.nii.gz ' outputfolder '/stats/Subject -t']);

files = dir([outputfolder 'stats/Subject*.nii.gz']);

for subject = 1:length(subjects)
    
    subjid = subjects_ordered{subject};
    
    eval(['!fslchfiletype NIFTI ' outputfolder 'stats/' files(subject).name ' ' outputfolder 'stats/' subjid '_skeletonized']);
    
end

delete([outputfolder 'stats/Subject*.nii.gz']);
delete([outputfolder 'stats/indsubs_FA_skeletonized.nii.gz']);

