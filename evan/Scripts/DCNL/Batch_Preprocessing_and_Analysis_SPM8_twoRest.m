%Batch_Preprocessing_and_Analysis_SPM8.m
%
% Given a list of subjects and an SPM8 batch template file, this script
% will automatically create a batch file for every subject and run them all
% in SPM8.
% 
% Written by E. Gordon, 12/2010

%List of subjects to run
subs = {'420'};
    %{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','283','281','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417' 


warning off

template_file = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/SPM8_Preprocessing_Template_twoRest.mat';

%subject loop
for subnum = 1:length(subs)
    

    
    subj = subs{subnum};
    
    disp(['Subject ' subj])
    
    cd(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/']);
    
    Nback_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/Nback/'];
    Rest_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/Rest/'];
    FirstRest_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/FirstRest/'];
    
    %automatically find raw images in the folder
    RestImages = dir([Rest_datafolder 'D*.img']);
    if isempty(RestImages)
        FourDimage = dir([Rest_datafolder '*.nii']);
        eval(['!fslsplit ' Rest_datafolder FourDimage.name ' ' Rest_datafolder 'DCNL-CV-' subj '- -t'])
        RestImages = dir([Rest_datafolder 'D*.nii.gz']);
        for i = 1:length(RestImages)
            if i<=2
                delete([Rest_datafolder RestImages(i).name]);
            else
                tounzip{i-2} = [Rest_datafolder RestImages(i).name];
            end
        end
        gunzip(tounzip)
        RestImages = dir([Rest_datafolder 'D*.nii']);
        delete([Rest_datafolder '*.nii.gz']);
    end
    clear tounzip
    
    FirstRestImages = dir([FirstRest_datafolder 'D*.img']);
    if isempty(FirstRestImages)
        FourDimage = dir([FirstRest_datafolder '*.nii']);
        eval(['!fslsplit ' FirstRest_datafolder FourDimage.name ' ' FirstRest_datafolder 'DCNL-CV-' subj '- -t'])
        FirstRestImages = dir([FirstRest_datafolder 'D*.nii.gz']);
        for i = 1:length(FirstRestImages)
            if i<=2
                delete([FirstRest_datafolder FirstRestImages(i).name]);
            else
                tounzip{i-2} = [FirstRest_datafolder FirstRestImages(i).name];
            end
        end
        gunzip(tounzip)
        FirstRestImages = dir([FirstRest_datafolder 'D*.nii']);
        delete([FirstRest_datafolder '*.nii.gz']);
    end
    clear tounzip
    
        
    NbackImages = dir([Nback_datafolder 'D*.img']);
    if isempty(NbackImages)
        FourDimage = dir([Nback_datafolder '*.nii']);
        eval(['!fslsplit ' Nback_datafolder FourDimage.name ' ' Nback_datafolder 'DCNL-CV-' subj '- -t'])
        NbackImages = dir([Nback_datafolder 'D*.nii.gz']);
        for i = 1:length(NbackImages)
            if i<=2
                delete([Nback_datafolder NbackImages(i).name]);
            else
                tounzip{i-2} = [Nback_datafolder NbackImages(i).name];
            end
        end
        gunzip(tounzip)
        NbackImages = dir([Nback_datafolder 'D*.nii']);
        delete([Nback_datafolder '*.nii.gz']);
    end
    clear tounzip
    
    
    Nback_statsoutputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/'];
    mkdir(Nback_statsoutputfolder);
    
    
    clear matlabbatch
    load(template_file);
    
    
    for scan = 1:length(NbackImages)
        matlabbatch{1}.spm.spatial.realign.estimate.data{1}{scan} = [Nback_datafolder NbackImages(scan).name ',1'];
    end
    matlabbatch{3}.spm.spatial.normalise.estwrite.subj.source{1} = [Nback_datafolder 'a' NbackImages(1).name ',1'];
    matlabbatch{5}.spm.stats.fmri_spec.dir{1} = Nback_statsoutputfolder;
    
    for scan = 1:length(RestImages)
        matlabbatch{8}.spm.spatial.realign.estimate.data{1}{scan} = [Rest_datafolder RestImages(scan).name ',1'];
    end
    matlabbatch{10}.spm.spatial.normalise.estwrite.subj.source{1} = [Rest_datafolder 'a' RestImages(1).name ',1'];
    
    for scan = 1:length(FirstRestImages)
        matlabbatch{12}.spm.spatial.realign.estimate.data{1}{scan} = [FirstRest_datafolder FirstRestImages(scan).name ',1'];
    end
    matlabbatch{14}.spm.spatial.normalise.estwrite.subj.source{1} = [FirstRest_datafolder 'a' FirstRestImages(1).name ',1'];
    
        
    batchfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8_Preprocessing.mat'];

    save(batchfilename, 'matlabbatch');

    spm_jobman('run',batchfilename)
    
end