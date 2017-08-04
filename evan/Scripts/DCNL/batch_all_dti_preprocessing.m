%batch_all_dti_preprocessing
%
%For any number of subjects, conducts all necessary preprocessing of
%Dicommed DTI data, including merging runs, eddy correction, DTIfit,
%bedposting, and calculation of registration parameters to standard space.
%
%Subjects/paths are specified at the top of the script
%
%You must have FSL loaded on your computer for this to run.
%
%Note that this takes ~7 hours per subject to run.
%
%Created by E. Gordon 10/09


%USER INPUT (and more down below)
%--------------------------------------------------------------------------

%Names of subjects
subjects = {'309'};
    %'343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
    %'101','102','110','113','118','120','122','125','127','132'};
    %,'138','147','150','151','154','156','159','160','161','162','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%'101',

%Run names (i.e. names of raw data folders).  These runs will be merged for
%each subject, if they're both ok. NOTE: at present this script can only
%handle merging of two runs. 
runs = {'run1' 'run2'};

%Runs to be used, for each subject.  This is a cell array in the same
%order as the subject names above, with each cell containing 'both', '1',
%or '2' to indicate which of the runs specified above is ok for that
%subject to use.  If both runs are not ok, no merging will take place, and
%only one run will be processed.
runstouse = {'both','both','both','both','both','both','both','both','both','both','both'};

%Base name of output file (to be written into the same folder the original
%data folders are in)
outputname = 'bothruns';

%template 2-run bvals and bvecs files
bvals2 = '/fmri/data3/Evan/Gene-Rest-Nback/Data/Temp/bvals_tworuns';
bvecs2 = '/fmri/data3/Evan/Gene-Rest-Nback/Data/Temp/bvecs_tworuns';

%template 1-run bvals and bvecs files
bvals1 = '/fmri/data3/Evan/Gene-Rest-Nback/Data/Temp/bvals';
bvecs1 = '/fmri/data3/Evan/Gene-Rest-Nback/Data/Temp/bvecs';

%define the standard brain to use
st_brain = '/data/apps/fsl/4.1.7/data/standard/MNI152_T1_2mm_brain';

%END USER INPUT
%--------------------------------------------------------------------------


for subnum = 1:length(subjects)
    subject = subjects{subnum};
    
    %USER INPUT
    %--------------------------------------------------------------------------
    
    %Location of DTI data folders specified at the top
    directory = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subject '/DTI/'];
    
    %END USER INPUT
    %--------------------------------------------------------------------------
    
    %Overwrite previous results
    try rmdir([directory outputname],'s');catch;end
    try rmdir([directory outputname '.bedpostX'],'s');catch;end
    
    %Make output directory
    mkdir([directory outputname]);
    
    %CD to that directory
    cd([directory outputname]);
    
%     %Check if we're merging data
%     if strcmp(runstouse{subnum},'both')
        
        %Concatenate the two runs
        %eval(['!fslmerge -t ' outputname ' ' directory runs{1} '/*' subject '*.nii ' directory runs{2} '/*' subject '*.nii']);
        eval(['!fslmerge -t ' outputname ' ' directory runs{1} '/*ep2d*.nii ' directory runs{2} '/*ep2d*.nii']);
        
        %Copy the 2-run bval and bvec template files
        copyfile(bvals2, [directory outputname '/bvals']);
        copyfile(bvecs2, [directory outputname '/bvecs']);
        
%     else
%         
%         %Just copy the good run
%         copyfile([directory runs{num2str(runstouse{subnum})} '/*' subject '*.nii'],[directory outputname '/' outputname]);
%         
%         %Turn it into a zipped nifti
%         eval(['!fslchfiletype NIFTI_GZ ' directory outputname '/' outputname]);
%         
%         %Copy the 1-run bval and bvec files
%         copyfile(bvals1, [directory outputname '/bvals']);
%         copyfile(bvecs1, [directory outputname '/bvecs']);
%         
%     end
    
    %Define the merged (or unmerged) data file
    input = [outputname '.nii.gz'];
    
    %Check to be sure orientation is radiological
    eval(['!fslorient -copyqform2sform ' input]);
    eval(['!fslorient -getorient ' input]);
    
    %Save the b0 data
    copyfile(input, ['copy_' input]);
    eval(['!fslsplit copy_' input]);
    delete(['copy_' input]);
    copyfile('vol0000.nii.gz','nodif.nii.gz');
    delete('vol*');
    
    %Conduct brain extraction
    eval(['!bet ./nodif.nii.gz ./nodif_brain -f 0.25 -g 0.0 -m']);
    
    %Conduct eddy correction
    eval(['!eddy_correct ./' input ' ./data 0']);
    
    %Conduct DTIfit
    eval(['!dtifit --data=./data.nii.gz --out=./dti --mask=./nodif_brain_mask.nii.gz --bvecs=./bvecs --bvals=./bvals'])
    
    %Conduct bedposting
    cd(directory)
    eval(['!bedpostx ' directory outputname ' -n 2']);
    
    %Register to standard space, but don't change data--just write results
    %to a transformation matrix
    eval(['!flirt -in ' directory outputname '.bedpostX/nodif_brain -ref ' st_brain ' -omat ' directory outputname '.bedpostX/xfms/diff2standard.mat -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12']);
    eval(['!convert_xfm -omat ' directory outputname '.bedpostX/xfms/standard2diff.mat -inverse ' directory outputname '.bedpostX/xfms/diff2standard.mat']);
    
    
    
end


