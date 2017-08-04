%Batch_VBM_allsteps
%
%For any number of subjects, registers all subjects to a common space,
%conducts VBM analysis within DARTEL, and thresholds the outputs.
%
%Subjects/paths are specified at the top of the script
%
%You must have FSL loaded on your computer, and SPM8 in your path for this to run
%
%
%Created by E. Gordon 10/09


%USER INPUT (and more down below)
%--------------------------------------------------------------------------

%Subject names
subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327'};

%Location of DARTEL batch template
DARTELtemplate = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/VBM_BATCH.mat';

%Location of registration template
Registertemplate = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/VBM_BATCH_register.mat';

%Location for DARTEL group analysis to be conducted (really this is just a
%place where the batch file can be saved)
DARTELlocation = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/VBM/';

%Threshold to use (.3 is apparently kind of standard?)
threshold = .3;

%END USER INPUT
%--------------------------------------------------------------------------


%Overwrite previous DARTEL batches
try rmfile(DARTELlocation,'s'); catch; end
mkdir(DARTELlocation);

%Load the DARTEL batch template
load(DARTELtemplate);

%Save the batch information
matlabbatch_DARTEL = matlabbatch;

clear matlabbatch;

currentdir = pwd;


for subnum = 1:length(subjects)
    
    subject = subjects{subnum};
    disp(['Registering ' subject])
    
    
    %USER INPUT (and more down below)
    %--------------------------------------------------------------------------
    
    %Location of subject's raw MPRAGE
    directory{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subject '/Struct/'];
    
    %Pattern that can be used to find this subject's MPRAGE
    MPRAGEpattern = 'DCNL*MPRAGE*';
    
    %END USER INPUT
    %--------------------------------------------------------------------------
    
    %Find the raw MPRAGE based on the pattern above
    orig_mprage = dir([directory MPRAGEpattern]);
    
    %Copy that image and header so as not to alter it
    copyfile([directory orig_mprage(1).name],[directory 'Registered_MPRAGE.hdr']);
    copyfile([directory orig_mprage(2).name],[directory 'Registered_MPRAGE.img']);
    
    %Load the registration batch template
    load(Registertemplate);
    
    %Put the found data into the registration batch
    matlabbatch{1}.spm.spatial.coreg.estimate.source{1} = [directory 'Registered_MPRAGE.img,1'];
    
    %Save and run the registration batch
    save([directory 'VBM_register.mat'], 'matlabbatch');
    spm_jobman('run',[directory 'VBM_register.mat']);
    
    clear matlabbatch
    
    %Put the results into the DARTEL template
    matlabbatch_DARTEL{1}.spm.tools.preproc8.channel.vols{subnum,1} = [directory 'Registered_MPRAGE.img,1'];
    
end

disp('Running DARTEL')
matlabbatch = matlabbatch_DARTEL;

%Save and run the DARTEL batch
save([DARTELlocation 'DARTEL.mat'], 'matlabbatch');
cd(DARTELlocation)
spm_jobman('run','DARTEL.mat');


cd(currentdir)

%Threshold the DARTEL outputs
disp('Thresholding outputs');

for subnum = 1:length(subjects)
    
    subject = subjects{subnum};
    
    %Threshold the GM ouput and turn it back into a Nifti
    eval(['!fslmaths ' directory{subnum} 'smwrc1Registered_MPRAGE -thr ' threshold ' ' directory{subnum} 'smwrc1Registered_MPRAGE_thr']);
    eval(['!fslchfiletype NIFTI ' directory{subnum} 'smwrc1Registered_MPRAGE_thr']);
    
    %Threshold the WM ouput and turn it back into a Nifti
    eval(['!fslmaths ' directory{subnum} 'smwrc2Registered_MPRAGE -thr ' threshold ' ' directory{subnum} 'smwrc2Registered_MPRAGE_thr']);
    eval(['!fslchfiletype NIFTI ' directory{subnum} 'smwrc2Registered_MPRAGE_thr']);
    
end