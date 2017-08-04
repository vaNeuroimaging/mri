%% Parameters

template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz';

T1string = 'T1';
T2string = 'T2';

normalized_voxdim = 3;

functional_sequences = {'restingstate'};

functional_sequence_TRs = {2.2};

functional_sequence_rawsize = {[64 64 34]};

fcprocess_sequences = {true};

skipframes = 5;




%% Steps to run

copy_rawdata = false;
T1_preprocess = false;
T1_segment = false;
T2_preprocess = false;
surface_registration = false;
BOLD_preprocess = true;
BOLD_fcprocess = false;
cifti_creation = false;


%% Get subjects


subjects = '/home/data/subjects/processing_list_111615.txt';
%subjects = {'MAV014'};

if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end


%% Set up parallel processing
warning off

maxworkers = 20;

subcount = length(subjects);
nworkers = min([subcount maxworkers]);

seq_counter = cell(subcount,1);


processingpool = parpool(nworkers);

parfor subnum = 1:subcount
%for subnum = 1:subcount
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    
    
    try
        
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
        
        if copy_rawdata
            disp(['Subject ' subject ': finding and copying data']);
            dlmwrite(logfile,'copying raw data to preprocessing folder','-append','delimiter','')
            
            T1_concatstr = ['fslmerge -t ' preprocfolder '/T1_concat'];
            T2_concatstr = ['fslmerge -t ' preprocfolder '/T2_concat'];
            
            firstT1 = true;
            
            for sess = 1:length(sessions)
                
                T1files = dir([rawfolder '/' sessions(sess).name '/o*' T1string '*']);
                for file = 1:length(T1files)
                    filename = [rawfolder '/' sessions(sess).name '/' T1files(file).name];
                    [~,sizestr] = system(['fslsize ' filename ' -s']);
                    tokens = tokenize(sizestr,' ');
                    datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                    if firstT1
                        T1size = datasize;
                        T1_concatstr = [T1_concatstr ' ' filename];
                        firstT1 = false;
                    else
                        if all(T1size==datasize)
                            T1_concatstr = [T1_concatstr ' ' filename];
                        end
                    end
                end
                
                T2files = dir([rawfolder '/' sessions(sess).name '/o*' T2string '*']);
                for file = 1:length(T2files)
                    filename = [rawfolder '/' sessions(sess).name '/' T2files(file).name];
                    T2_concatstr = [T2_concatstr ' ' filename];
                end
                
                %Find all BOLD runs specified
                
                                                                
                for seq = 1:length(functional_sequences)
                    
                    thisseq_files = dir([rawfolder '/' sessions(sess).name '/*' functional_sequences{seq} '*.nii.gz']);
                    for file = 1:length(thisseq_files)
                        filename = [rawfolder '/' sessions(sess).name '/' thisseq_files(file).name];
                        
                        %check if data is the right size
                        [~,sizestr] = system(['fslsize ' filename ' -s']);
                         tokens = tokenize(sizestr,' ');
                         datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                         
                         if all(datasize(1:3)==functional_sequence_rawsize{seq}) && datasize(4) > 20
                        
                            %load data, trim first few frames
                            data = load_untouch_nii(filename);
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
                
                
                
                
                
                 
%                 %Find all T1s and prepare to align and average them
%                 T1folders = dir([rawfolder '/' sessions(sess).name '/NIFTI/*T1*']);
%                 for folder = 1:length(T1folders);
%                     filename = dir([rawfolder '/' sessions(sess).name '/NIFTI/' T1folders(folder).name '/*.nii.gz']);
%                     T1_concatstr = [T1_concatstr ' ' rawfolder '/' sessions(sess).name '/NIFTI/' T1folders(folder).name '/' filename(1).name];
%                 end
%                 
%                 %Find all T2s and prepare to align and average them
%                 T2folders = dir([rawfolder '/' sessions(sess).name '/NIFTI/*T2*']);
%                 for folder = 1:length(T2folders);
%                     filename = dir([rawfolder '/' sessions(sess).name '/NIFTI/' T2folders(folder).name '/*.nii.gz']);
%                     T2_concatstr = [T2_concatstr ' ' rawfolder '/' sessions(sess).name '/NIFTI/' T2folders(folder).name '/' filename(1).name];
%                 end
%                 
%                 %Find all BOLD runs specified
%                 disp(['Subject ' subject ': finding and copying BOLD runs']);
%                                                                 
%                 for seq = 1:length(functional_sequences)
%                     thisseq_folders = dir([rawfolder '/' sessions(sess).name '/NIFTI/*' functional_sequences{seq} '*']);
%                     for folder = 1:length(thisseq_folders);
%                         filename = dir([rawfolder '/' sessions(sess).name '/NIFTI/' thisseq_folders(folder).name '/*.nii.gz']);
%                         
%                         %keep track of session names, session number for each sequence, and overall session number
%                         seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
%                         BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
%                         sesscounter(end+1) = sess;
%                         
%                         %load data, trim first few frames, and save to preprocessing folder
%                         data = load_untouch_nii([rawfolder '/' sessions(sess).name '/NIFTI/' thisseq_folders(folder).name '/' filename(1).name]);
%                         data.img = data.img(:,:,:,skipframes+1:end);
%                         %data.img = data.img(:,:,:,skipframes+1:end-12);
%                         data.hdr.dime.dim(5) = size(data.img,4);
%                         save_untouch_nii(data,[BOLDrunnames{end} '.nii.gz'])
%                         
%                     end
%                 end
                
            end
            
            %T1s: concatenate, align, and average
            disp(['Subject ' subject ': aligning and averaging both T1s and T2s']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(T1_concatstr);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' preprocfolder '/T1_concat -refvol 0']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' preprocfolder '/T1_concat_mcf -Tmean ' preprocfolder 'T1_avg']);
            delete([preprocfolder '/T1_concat.nii.gz'])
            delete([preprocfolder '/T1_concat_mcf.nii.gz'])
            
            %T2s: concatenate, align, and average
            [systemresult{end+1,1},systemresult{end+1,2}] = system(T2_concatstr);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' preprocfolder '/T2_concat -refvol 0']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' preprocfolder '/T2_concat_mcf -Tmean ' preprocfolder 'T2_avg']);
            delete([preprocfolder '/T2_concat.nii.gz'])
            delete([preprocfolder '/T2_concat_mcf.nii.gz'])
            
        else
            %keep track of functional session names, session number for each sequence, and overall session number
            for sess = 1:length(sessions)
                for seq = 1:length(functional_sequences)
                     thisseq_files = dir([rawfolder '/' sessions(sess).name '/*' functional_sequences{seq} '*.nii.gz']);
                     for file = 1:length(thisseq_files)
                         filename = [rawfolder '/' sessions(sess).name '/' thisseq_files(file).name];
                         
                         [~,sizestr] = system(['fslsize ' filename ' -s']);
                         tokens = tokenize(sizestr,' ');
                         datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                         
                         if all(datasize(1:3)==functional_sequence_rawsize{seq}) && datasize(4) > 20
                             %keep track of session names, session number for each sequence, and overall session number
                             seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
                             BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
                             sesscounter(end+1) = sess;
                        end
                    end
                end
            end
            
            
        end
        
        
        
        
        %% T1 Preprocessing
        cd(preprocfolder)
        
        T1file = 'T1_avg';
        
        if T1_preprocess
            
            dlmwrite(logfile,'preprocessing T1','-append','delimiter','')
            
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting, and skull-stripping T1, and registering to MNI']);
            
            %reorient, crop, and bias-correct T1 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T1file ' -o struct --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            copyfile([preprocfolder '/struct.anat/T1_biascorr.nii.gz'],[preprocfolder '/' T1file '_biascorr.nii.gz'])
             T1file = [T1file '_biascorr'];
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['chmod 777 -R ' preprocfolder '/struct.anat'])
            delete([preprocfolder '/struct.anat/*']);
            rmdir([preprocfolder '/struct.anat']);
            
            %brain-extract T1 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' T1file ' ' T1file '_bet -R -g -.1']);
            
            %register T1 to MNI
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file '_bet -interp spline -dof 12 -ref ' template ' -omat T1_2MNI.mat -out ' T1file '_bet_MNI']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file ' -interp spline -ref ' template ' -applyxfm -init T1_2MNI.mat -out ' T1file '_MNI']);
            
        else
            
            %track name
            T1file = [T1file '_biascorr'];
            
        end
        
        
        %% T1 Segmentation
        
        if T1_segment
            
            disp(['Subject ' subject ': segmenting T1']);
            dlmwrite(logfile,'segmenting T1','-append','delimiter','')
            
            try [~,~,~] = rmdir(freesurferfolder,'s'); catch; end
            
            mkdir(freesurferfolder)
            
            %copy bias-corrected T1 to freesurfer folder
            copyfile([preprocfolder '/' T1file '.nii.gz'],[freesurferfolder '/' T1file '.nii.gz'])
            cd(freesurferfolder)
            
            %transform T1 to .mgz format
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['recon-all -subjid ' subject ' -sd ' freesurferfolder ' -i ' freesurferfolder '/' T1file '.nii.gz']);
            %segment T1
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['recon-all -subjid ' subject ' -sd ' freesurferfolder ' -all']);
            
            %move files around
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' freesurferfolder '/' subject '/* ' freesurferfolder]);
            rmdir([freesurferfolder '/' subject ])
            cd(preprocfolder)
            
            if ~systemresult{end,1}
                %get segmentation masks into nifti format
                systemresult = make_fs_masks(freesurferfolder,normalized_voxdim,[preprocfolder '/T1_2MNI.mat'],template,systemresult);
                copyfile([freesurferfolder '/nusmask/aparc+aseg_cerebralwm.nii.gz'],preprocfolder)
            else
                error('Problem detected with Freesurfer')
            end
            
        end
        
        
        %% T2 Preprocessing
        
        cd(preprocfolder)
        
        T2file = 'T2_avg';
        
        if T2_preprocess
            
            dlmwrite(logfile,'preprocessing T2','-append','delimiter','')
            
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting, and skull stripping T2, and registering to T1 and through to MNI']);
            
            %reorient, crop, and bias-correct T2 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T2file ' -o struct -t T2 --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            copyfile([preprocfolder '/struct.anat/T2_biascorr.nii.gz'],[preprocfolder '/' T2file '_biascorr.nii.gz'])
            T2file = [T2file '_biascorr'];
            
            delete([preprocfolder '/struct.anat/*']);
            rmdir([preprocfolder '/struct.anat']);
            
            %brain-extract T2 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' T2file ' ' T2file '_bet -R -g -.1']);
            
            %register T2 to T1, and concatenate with T1 to MNI
            %OPTION 1
            %[status(end+1) result{end+1}] = system(['flirt -in ' T2file ' -ref ' T1file ' -interp spline -dof 6 -omat T2_2T1.mat -out T2_T1space']);
            
            
            %OPTION 2
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -dof 6 -omat T2_2T1_init.mat']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -interp spline -dof 6 -cost bbr -wmseg aparc+aseg_cerebralwm -init T2_2T1_init.mat -omat T2_2T1.mat -out T2_T1space -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
            
            [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat T1_2T2.mat -inverse T2_2T1.mat');
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in aparc+aseg_cerebralwm -ref ' T2file ' -applyxfm -init T1_2T2.mat -interp nearestneighbour -out aparc+aseg_cerebralwm_T2space']);
            
            [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat T2_2MNI.mat -concat T1_2MNI.mat T2_2T1.mat');
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' template ' -applyxfm -init T2_2MNI.mat -out ' T2file '_MNI']);
            
        else
            
            %track T2 name
            T2file = [T2file '_biascorr'];
            
        end
        
        
        %% Surface-based registration
        
        if surface_registration
            
            disp(['Subject ' subject ': surface-based registration']);
            dlmwrite(logfile,'conducting surface registration','-append','delimiter','')
            
            mkdir(fsLRfolder)
            
            prevsize = size(systemresult,1);
            systemresult = PostFreeSurferPipeline_fsavg2fslr_long_singlesub(subject,T1file,freesurferfolder,fsLRfolder,preprocfolder,systemresult);
            if any(cell2mat(systemresult(prevsize+1:end,1)))
                error('Problem detected with surface registration')
            end
            
        end
        
        
        %% Functional preprocessing
        
        cd(preprocfolder)
        
        if BOLD_preprocess
            
            dlmwrite(logfile,'preprocessing BOLD data','-append','delimiter','')
            
            %make a file that tracks session number
            sess_tracker = [preprocfolder '/BOLDruns_sessions.txt'];
            delete(sess_tracker);
            fid = fopen(sess_tracker,'at'); %open the output file for writing
            fclose(fid);
            
            brainmask = load_untouch_nii([freesurferfolder '/nusmask/aparc+aseg_brainmask_mask_' num2str(normalized_voxdim) num2str(normalized_voxdim) num2str(normalized_voxdim) '.nii.gz']);
            
            mergestr = 'fslmerge -t allBOLDavgs ';
            
            disp(['Subject ' subject ': slice time correction and motion correction']);
            
            ign = rmdir('allBOLDavgs_mcf.mat*');
            
            for runnum = 1:length(BOLDrunnames)
                BOLDfile = BOLDrunnames{runnum};
                
                %get short filename
                [~,BOLDfilename,~] = fileparts(BOLDfile);
                
                %get TR for this run
                for i = 1:length(functional_sequences)
                    if strcmp(BOLDfilename(1:length(functional_sequences{i})),functional_sequences{i});
                        TR = functional_sequence_TRs{i};
                    end
                end
                
                %slice-time correct BOLD images
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['slicetimer -i ' BOLDfile ' -o ' BOLDfile '_st -r ' num2str(TR) ' --ocustom=/home/data/scripts/Resources/sliceorder.txt']);
                
                ign = rmdir([BOLDfile '_st_mcf.mat*']);
                %motion-correct BOLD images to first image to get motion plots
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' BOLDfile '_st -refvol 0 -plots -mats']);
                
                %split slice-time corrected timeseries into separate volumes
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslsplit ' BOLDfile '_st ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf -t']);
                
                %get avg BOLD for this run
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' BOLDfile '_st_mcf -Tmean ' BOLDfile '_st_mcf_avg']);
                
                %Add to string to concatenate avg BOLD images
                mergestr = [mergestr ' ' BOLDfile '_st_mcf_avg'];
                
                %delete motion corrected image--we will include this in a single-step atlas transformation below
                delete([BOLDfile '_st_mcf.nii.gz'])
                
            end
            
            %align average BOLD volumes across runsm calculate overall average BOLD volume, and brain-extract
            [systemresult{end+1,1},systemresult{end+1,2}] = system(mergestr);
            [systemresult{end+1,1},systemresult{end+1,2}] = system('mcflirt -in allBOLDavgs -dof 9 -mats');
            [systemresult{end+1,1},systemresult{end+1,2}] = system('fslmaths allBOLDavgs_mcf -Tmean BOLD_avg');
            [systemresult{end+1,1},systemresult{end+1,2}] = system('bet BOLD_avg BOLD_avg_bet -R');
            
            
            %register overall average BOLD to T1 and concatenate through to MNI
            disp(['Subject ' subject ': register average BOLD to T1 and through to MNI, apply to individual runs, and normalize intensity to mode 1000']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in BOLD_avg_bet -ref '  T1file '_bet -dof 12 -omat BOLD_avg_2T1_init.mat']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in BOLD_avg_bet -ref '  T1file '_bet -dof 12 -cost bbr -wmseg aparc+aseg_cerebralwm -init BOLD_avg_2T1_init.mat -omat BOLD_avg_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat BOLD_avg_2MNI.mat -concat T1_2MNI.mat BOLD_avg_2T1.mat']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in BOLD_avg -interp spline -ref ' template ' -applyisoxfm ' num2str(normalized_voxdim) ' -init BOLD_avg_2MNI.mat -out BOLD_avg_MNI']);
            
            avgBOLDmats = dir('allBOLDavgs_mcf.mat/MAT*');
            
            for runnum = 1:length(BOLDrunnames)
                BOLDfile = BOLDrunnames{runnum};
                
                %get short filename
                [~,BOLDfilename,~] = fileparts(BOLDfile);
                
                %get this run to avg run transform
                BOLD_to_overallBOLD_mat = ['allBOLDavgs_mcf.mat/' avgBOLDmats(runnum).name];
                
                %concatenate motion correction transform with reference volume to MNI transform
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_avg_2MNI.mat -concat BOLD_avg_2MNI.mat ' BOLD_to_overallBOLD_mat]);
                
                %get volume numbers
                vols = dir([BOLDfile '_st_mcf.mat/MAT*']);
                fslmergestr = ['fslmerge -t ' BOLDfile '_st_mcf_MNI'];
                %for each volume
                for i = 1:length(vols)
                    %concatenate motion correction transform with run avg to MNI transform
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat -concat  ' BOLDfile '_avg_2MNI.mat ' BOLDfile '_st_mcf.mat/' vols(i).name]);
                    
                    %apply concatenated transform
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf' vols(i).name(end-3:end) ' -interp spline -ref ' template ' -applyisoxfm ' num2str(normalized_voxdim) ' -init '  BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat -out ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf_MNI' vols(i).name(end-3:end)]);
                    
                    %prepare to merge all volumes back into a 4D image
                    fslmergestr = [fslmergestr ' ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf_MNI' vols(i).name(end-3:end)];
                end
                %merge volumes
                [systemresult{end+1,1},systemresult{end+1,2}] = system(fslmergestr);
                
                %remove unneeded files
                delete([BOLDfile '_st_mcf.mat/*'])
                rmdir([BOLDfile '_st_mcf.mat'])
                delete([BOLDfile '_st.nii.gz'])
                delete([BOLDfile '_2T1_init.mat'])
                
                
                %Intensity normalization to mode within-brain value of 1000
                data = load_untouch_nii([BOLDfile '_st_mcf_MNI.nii.gz']);
                data_inmask = data.img .* repmat(brainmask.img,[1,1,1,size(data.img,4)]) .* (data.img > 100);
                data_inmask = data_inmask(data_inmask>0);
                [counts,edges] = histcounts(data_inmask,1000);
                [~,maxind] = max(counts);
                modeval = mean([edges(maxind) edges(maxind+1)]);
                data.img = data.img + (1000 - modeval);
                save_untouch_nii(data,[BOLDfile '_st_mcf_MNI.nii.gz']);
                data = []; data_inmask = [];
                %[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' BOLDfile '_st_mcf_MNI -thr 100 -ing 1000 ' BOLDfile '_st_mcf_MNI']);
                
                %track session number
                dlmwrite(sess_tracker,[BOLDfile '_st_mcf_MNI.nii.gz ' num2str(sesscounter(runnum))],'-append','delimiter','');%write the data to the output file
                
            end
            
            delete('allBOLDavgs_mcf.mat/*')
            rmdir('allBOLDavgs_mcf.mat')
            brainmask = [];
            
        end
        
        
        %% Fc-processing
        
        
        if BOLD_fcprocess
            
            mkdir(fc_processfolder)
            
            dlmwrite(logfile,'fc-processing BOLD data','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                if fcprocess_sequences{seq}
                    
                    disp(['Subject ' subject ', ' functional_sequences{seq} ': fc-processing']);
                    %FCPROCESS
                    
                end
            end
            
        end
        
        %% Cifti creation
        
        if cifti_creation
            
            mkdir(ciftifolder)
            
            dlmwrite(logfile,'creating ciftis','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': mapping functional data to cortical surface']);
                
                if fcprocess_sequences{seq}
                    
                    systemresult = post_fc_processing_batch_singlesub(subject,voldata_tomap,ciftifolder,fsLRfolder,functional_sequence_TRs{seq},functional_sequences{seq},tmaskfile,preproc_datalist,'/home/data/scripts/Cifti_creation/post_fc_processing_batch_params_singlesub.m',systemresult);
                    
                else
                    
                    BOLDmergedfile = [preprocfolder '/' functional_sequences{seq} '_merged'];
                    
                    fslmergestr = ['fslmerge -t ' BOLDmergedfile ' '];
                    for runnum = 1:length(BOLDrunnames)
                        [~,BOLDfilename,~] = fileparts(BOLDrunnames{runnum});
                        if strcmp(BOLDfilename(1:length(functional_sequences{seq})),functional_sequences{seq});
                            fslmergestr = [fslmergestr BOLDrunnames{runnum} '_st_mcf_MNI.nii.gz '];
                        end
                    end
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(fslmergestr);
                    
                    systemresult = post_fc_processing_batch_singlesub(subject,BOLDmergedfile,ciftifolder,fsLRfolder,functional_sequence_TRs{seq},functional_sequences{seq},[],[],'/home/data/scripts/Cifti_creation/post_fc_processing_batch_params_singlesub.m',systemresult);
                    
                    delete(BOLDmergedfile)
                    
                end
            end
            
        end
        
        
        disp(['Subject ' subject ' COMPLETED all processing steps without error!'])
        dlmwrite(logfile,'all processing complete.','-append','delimiter','')
        
    catch errorinfo
        
        dlmwrite(logfile,'ERROR:','-append','delimiter','')
        
        errormessage = [errorinfo.stack.file ', line ' num2str(errorinfo.stack.line) ': ' errorinfo.message];
        
        probleminds = find(cell2mat(systemresult(:,1)));
        for i = 1:length(probleminds)
            dlmwrite(logfile,systemresult{probleminds(i),2},'-append','delimiter','');
        end
        dlmwrite(logfile,errormessage,'-append','delimiter','');
        
        disp(['Subject ' subject ' FAILED to complete processing! Check error log ' logfile])
        
    end
    
    
end

delete(processingpool)