%Batch_Preprocessing_and_Analysis_SPM8.m
%
% Given a list of subjects and an SPM8 batch template file, this script
% will automatically create a batch file for every subject and run them all
% in SPM8.
% 
% Written by E. Gordon, 12/2010

%List of subjects to run
subs = {'281','292','301','307','322'};
    %101,'102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274'


warning off

template_file = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/SPM8_Preprocessing_Template.mat';

%subject loop
for subnum = 1:length(subs)
    

    
    subj = subs{subnum};
    
    disp(['Subject ' subj])
    
    cd(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/']);
    
    %define folder for this run
    orig_Nback_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Nback/'];
    Nback_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/Nback/'];
    
    mkdir([orig_Nback_datafolder '../SPM8/']);
    mkdir(Nback_datafolder);
    copyfile([orig_Nback_datafolder 'D*'],Nback_datafolder);
    
    orig_Rest_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Rest/'];
    Rest_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/Rest/'];
    
    mkdir([orig_Rest_datafolder '../SPM8/']);
    mkdir(Rest_datafolder);
    copyfile([orig_Rest_datafolder 'D*'],Rest_datafolder);
    
    
    
    
    %automatically find raw images in the folder
    NbackImages = dir([Nback_datafolder 'D*.img']);
    RestImages = dir([Rest_datafolder 'D*.img']);
    
    Nback_statsoutputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/'];
    mkdir(Nback_statsoutputfolder);
    
    
    clear matlabbatch
    load(template_file);
    
    
    for scan = 1:length(NbackImages)
        matlabbatch{1}.spm.spatial.realign.estimate.data{1}{scan} = [Nback_datafolder NbackImages(scan).name ',1'];
    end
    matlabbatch{3}.spm.spatial.normalise.estwrite.subj.source{1} = [Nback_datafolder 'a' NbackImages(1).name ',1'];
    matlabbatch{5}.spm.stats.fmri_spec.dir{1} = Nback_statsoutputfolder;
    
    if length(RestImages) > 0
        
        
        for scan = 1:length(RestImages)
            matlabbatch{8}.spm.spatial.realign.estimate.data{1}{scan} = [Rest_datafolder RestImages(scan).name ',1'];
        end
        matlabbatch{10}.spm.spatial.normalise.estwrite.subj.source{1} = [Rest_datafolder 'a' RestImages(1).name ',1'];
        
    else
        matlabbatch = matlabbatch(1:7);
    end
    
    
    
    batchfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8_Preprocessing.mat'];

    save(batchfilename, 'matlabbatch');

    spm_jobman('run',batchfilename)

    
end