subjects = {'MAV004'};
%{'MAV008','MAV009','MAV015','MAV016','MAV019','MAV020','MAV024','MAV029','MAV033','MAV036','MAV040','MAV042','MAV043','MAV044','MAV046','MAV047','MAV053','MAV055','MAV056','MAV058'};

functional_sequences = {'RSFC'};%,'Movie'};%

functional_sequence_TRs = {3.0};% 2.2

functional_sequence_rawsize = {[80 80 34]};%[64 64 34]

subcount = length(subjects);

seq_counter = cell(subcount,1);
QCtrackingfile = '/home/data/subjects/QC/MRI_scan_tracker.xlsx';
[~,~,QCdata]=xlsread(QCtrackingfile);


for subnum = 1:length(subjects)
    
    subject = subjects{subnum};
    
    rawfolder = ['/home/data/subjects/' subject '/raw/'];
    DTIfolder = ['/home/data/subjects/' subject '/DTI/'];
    preprocfolder = ['/home/data/subjects/' subject '/preprocessed/']; mkdir(preprocfolder)
    freesurferfolder = ['/home/data/subjects/' subject '/freesurfer/'];
    fsLRfolder = ['/home/data/subjects/' subject '/fs_LR/'];
    fc_processfolder = ['/home/data/subjects/' subject '/fc_processed/'];
    ciftifolder = ['/home/data/subjects/' subject '/cifti/'];
    logfolder = ['/home/data/subjects/' subject '/processing_logs/']; mkdir(logfolder)
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    sesscounter = [];
    
    sessions = dir([rawfolder '/']);
    folderindex = struct2cell(sessions); folderindex = logical(cell2mat(folderindex(4,:)));
    sessions = sessions(folderindex);
    sessions = sessions(3:end);
    sessnums_unsorted = zeros(1,length(sessions));
    for sess = 1:length(sessions)
        tokens = tokenize(sessions(sess).name,'-');
        sessnums_unsorted(sess) = str2num(tokens{end});
    end
    [sessnums, sessi] = sort(sessnums_unsorted,'ascend');
    sessions = sessions(sessi);
    
    
    %sessions = dir([rawfolder]); sessions = sessions(3);
    
    for seq = 1:length(functional_sequences)
        seq_counter{subnum}{seq} = 0;
    end
    BOLDrunnames = cell(0,1);
    
    subQCindices = find(strcmp(QCdata(:,1),subject));
    subQC_sessions = cell2mat(QCdata(subQCindices,2));
    
    runtrackerfile = [preprocfolder '/BOLD_run_tracker.txt'];
    delete(runtrackerfile);
    fid = fopen(runtrackerfile,'at');
    fclose(fid);
    
    
    
    for sess = 1:length(sessions)
        
        dashloc = strfind(sessions(sess).name,'-');
        sessnum = str2num(sessions(sess).name(dashloc+1:end));
        sessQCindices = subQCindices(subQC_sessions==sessnum);
        sessQC_decision = cell2mat(QCdata(sessQCindices,7));
        sessQC_acquisitionnum = cell2mat(QCdata(sessQCindices,8));
        
        %Find all valid BOLD runs
        for seq = 1:length(functional_sequences)
            
            thisseq_files = dir([rawfolder '/' sessions(sess).name '/*' functional_sequences{seq} '*.nii.gz']);
            for file = 1:length(thisseq_files)
                filename = [rawfolder '/' sessions(sess).name '/' thisseq_files(file).name];
                
                extpoint = strfind(filename,'.nii.gz');
                acquisitionnum = str2num(filename(extpoint-2:extpoint-1));
                
                %check if data is the right size
                [~,sizestr] = system(['fslsize ' filename ' -s']);
                tokens = tokenize(sizestr,' ');
                datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                correctsize = logical(all(datasize(1:3)==functional_sequence_rawsize{seq}));
                
                %check if data has many frames (i.e. wasn't aborted)
                data_islong = (datasize(4) > 20);
                
                %check if data has been flagged by manual QC check (or if it's not in the QC check file)
                if isempty(acquisitionnum)
                    flaggedasbad = false;
                else
                    flaggedasbad = (sessQC_decision(sessQC_acquisitionnum==acquisitionnum)==0);
                    if ~any(sessQC_acquisitionnum==acquisitionnum)
                        flaggedasbad = false;
                    end
                end
                
                if correctsize && data_islong && (~flaggedasbad)
                    
                    %keep track of session names, session number for each sequence, and overall session number
                    seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
                    BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
                    sesscounter(end+1) = sessnums(sess);
                    
                    string_towrite = [filename ' ' BOLDrunnames{end} '.nii.gz'];
                    dlmwrite(runtrackerfile,string_towrite,'-append','delimiter','')
                    
                end
                
            end
        end
        
    end
    
end