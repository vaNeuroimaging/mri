%Batch_Analysis_SPM8.m
%
% Given a list of subjects and an SPM8 batch template file, this script
% will automatically create a batch file for every subject and run them all
% in SPM8.
% 
% Written by E. Gordon, 12/2010

%List of subjects to run
subs = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
    %'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274',

%
%133
warning off

template_file = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/SPM8_Analysis_FirstRest.mat';

%subject loop
for subnum = 1:length(subs)
    

    
    subj = subs{subnum};
    
    disp(['Subject ' subj])
    
    Datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/FirstRest/'];
    
    Statsoutputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/FirstRest_forDCM/'];
    
    Images = dir([Datafolder 'fswa*.img']);
    if isempty(Images)
        Images = dir([Datafolder 'fsw*.nii']);
    end
    
    clear matlabbatch
    load(template_file);
    
    
    for scan = 1:length(Images)
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans{scan} = [Datafolder Images(scan).name ',1'];
    end
    matlabbatch{1}.spm.stats.fmri_spec.dir{1} = Statsoutputfolder;
    
    try rmdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'s'); catch; end
    mkdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1});
    
    batchfilename = [Statsoutputfolder(1:end-1) '.mat'];

    save(batchfilename, 'matlabbatch');

    spm_jobman('run',batchfilename)

    
end