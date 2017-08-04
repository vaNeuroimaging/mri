%% Parameters

template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz';

normalized_voxdim = 3;

functional_sequences = {'restingstate'};

functional_sequence_TRs = {2.2};

fcprocess_sequences = {false};

skipframes = 5;




%% Steps to run

copy_rawdata = false;
T1_preprocess = false;
T1_segment = true;
T2_preprocess = true;
surface_registration = true;
BOLD_preprocess = false;
BOLD_fcprocess = false;
cifti_creation = false;


%% Get subjects

%subjects = '/home/data/subjects/processing_list_102315.txt';
subjects = {'MAV013','MAV014','MAV020','MAV023','MAV025'};

if iscell(subjects)
    
elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end


%% Set up parallel processing

maxworkers = 24;

subcount = length(subjects);
nworkers = min([subcount maxworkers]);

seq_counter = cell(subcount,1);


for subnum = 1:subcount
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    warning off
    
    
    
    %% Set up folders
    
    
    rawfolder = ['/home/data/subjects/' subject '/raw/'];
    preprocfolder = ['/home/data/subjects/' subject '/preprocessed/']; mkdir(preprocfolder)
    freesurferfolder = ['/home/data/subjects/' subject '/freesurfer/'];
    fsLRfolder = ['/home/data/subjects/' subject '/fs_LR/'];
    fc_processfolder = ['/home/data/subjects/' subject '/fc_processed/'];
    ciftifolder = ['/home/data/subjects/' subject '/cifti/'];
    logfolder = ['/home/data/subjects/' subject '/processing_logs/']; mkdir(logfolder)
    
    logfile = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile,' '); logfile = [tokens{1} '_' tokens{2}];
    %Create error tracking file to record issues with processing
    fid = fopen(logfile,'at');
    fclose(fid);
    
    
    %% Find nifti data in raw folder and move to preprocessed folder
    
    
    sesscounter = [];
    
    sessions = dir([rawfolder '/']);
    folderindex = struct2cell(sessions); folderindex = logical(cell2mat(folderindex(4,:)));
    sessions = sessions(folderindex);
    sessions = sessions(3:end);
    %sessions = dir([rawfolder]); sessions = sessions(3);
    
    for seq = 1:length(functional_sequences)
        seq_counter{subnum}{seq} = 0;
    end
    BOLDrunnames = cell(0,1);
    
    
    disp(['Subject ' subject ': finding and copying data']);
    dlmwrite(logfile,'copying raw data to preprocessing folder','-append','delimiter','')
    
    T1_concatstr = ['fslmerge -t ' preprocfolder '/T1_concat'];
    T2_concatstr = ['fslmerge -t ' preprocfolder '/T2_concat'];
    
    for sess = 1:length(sessions)
        
        T1files = dir([rawfolder '/' sessions(sess).name '/o*T1*']);
        for file = 1:length(T1files)
            filename = [rawfolder '/' sessions(sess).name '/' T1files(file).name];
            T1_concatstr = [T1_concatstr ' ' filename];
        end
        
        T2files = dir([rawfolder '/' sessions(sess).name '/o*T2*']);
        for file = 1:length(T2files)
            filename = [rawfolder '/' sessions(sess).name '/' T2files(file).name];
            T2_concatstr = [T2_concatstr ' ' filename];
        end
        
        %Find all BOLD runs specified
        
        
        for seq = 1:length(functional_sequences)
            
            thisseq_files = dir([rawfolder '/' sessions(sess).name '/*' functional_sequences{seq} '*.nii.gz']);
            for file = 1:length(thisseq_files)
                filename = [rawfolder '/' sessions(sess).name '/' thisseq_files(file).name];
                
                
                
                
                %load data, check if it's the right size, and trim first few frames
                data = load_untouch_nii(filename);
                if all(data.hdr.dime.dim(2:4) == [64 64 34])
                    data.img = data.img(:,:,:,skipframes+1:end);
                    data.hdr.dime.dim(5) = size(data.img,4);
                    
                    %keep track of session names, session number for each sequence, and overall session number
                    seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
                    BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
                    sesscounter(end+1) = sess;
                    
                    %save to preprocessing folder
                    save_untouch_nii(data,[BOLDrunnames{end} '.nii.gz'])
                end
                
            end
        end
    end
end