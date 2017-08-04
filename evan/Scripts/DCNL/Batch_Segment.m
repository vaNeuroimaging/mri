%Batch_Preprocessing_and_Analysis.m
%
% Given a list of subjects and an SPM8 batch template file, this script
% will automatically create a batch file for every subject and run them all
% in SPM8.
% This is easiest to do for preprocessing, but it can be done for
% first-level statistical analysis as well (especially for analysis designs that do
% not vary from subject to subject).
%
% Written by E. Gordon, 04/2009

%List of subjects to run
subs = {'402'};
%{'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420','199','269','189','110'};

warning off

%subject loop
for subnum = 1:length(subs)
    subj = subs{subnum};
    disp(subj)
    folder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/'];
    %cd(folder);
    Imagesinfolder = dir([folder 'D*MPRAGE*.img']);
    if isempty(Imagesinfolder)
        Imagesinfolder = dir([folder 'S*MPRAGE*.nii']);
    end
    
    load /fmri/data3/Evan/Gene-Rest-Nback/Scripts/Segment_template.mat
    
    matlabbatch{1}.spm.spatial.preproc.data{1} = [folder Imagesinfolder(1).name ',1'];
    
    save(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Segment.mat'], 'matlabbatch');
 
%     %Call SPM and run it all
     spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Segment.mat']);
     
end

    
 