%Batch_Preprocessing_and_Analysis_SPM8.m
%
% Given a list of subjects and an SPM8 batch template file, this script
% will automatically create a batch file for every subject and run them all
% in SPM8.
% 
% Written by E. Gordon, 12/2010

%List of subjects to run
subs = {'181','182','110','189','199','269','101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};


warning off

template_file = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/SPM8_Nback_Analysis_only.mat';

%subject loop
for subnum = 1:length(subs)
    

    
    subj = subs{subnum};
    
    disp(['Subject ' subj])
    
    cd(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/']);
    
    %define folder for this run
    Nback_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/Nback/'];
    
    
    
    
    
    %automatically find raw images in the folder
    NbackImages = dir([Nback_datafolder 'sw*.img']);
    
    Nback_statsoutputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/'];
    try rmdir(Nback_statsoutputfolder,'s');catch;end
    mkdir(Nback_statsoutputfolder);
    
    
    clear matlabbatch
    load(template_file);
    
    
    for scan = 1:length(NbackImages)
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans{scan} = [Nback_datafolder NbackImages(scan).name ',1'];
    end
    matlabbatch{1}.spm.stats.fmri_spec.dir{1} = Nback_statsoutputfolder;
    
    try rmdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'s'); catch; end
    mkdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1});
    
    
    batchfilename = [Nback_statsoutputfolder(1:end-1) '.mat'];

    save(batchfilename, 'matlabbatch');

    spm_jobman('run',batchfilename)

    
end