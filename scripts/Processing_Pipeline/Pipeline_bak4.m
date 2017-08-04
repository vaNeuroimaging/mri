%% Subjects

%subjects = '/home/data/subjects/processing_list_10min.txt';
subjects = '/home/data/subjects/processing_list_010516.txt';
%subjects = '/home/data/subjects/list_forJDP.txt';
%subjects = {'MAV014'};



%% Steps to run

choose_T1 = true;
copy_rawBOLDdata = true;
T1_preprocess = false;
T1_segment = false;
T2_preprocess = false;
surface_registration = false;
BOLD_preprocess = false;
BOLD_fcprocess = true;
cifti_creation = true;
cifti_correlation = true;


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

FDthresh = .2;

mintime = 10; % minutes

sliceorderfile = '/home/data/scripts/Resources/sliceorder_interleavedPhillips.txt';

QCtrackingfile = '/home/data/subjects/QC/MRI_scan_tracker.xlsx';



%% Get subjects


if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end



subcount = length(subjects);

seq_counter = cell(subcount,1);

[~,~,QCdata]=xlsread(QCtrackingfile);

%% Choose best T1 image

if choose_T1
    for subnum = 1:subcount
        subject = subjects{subnum};
        
        rawfolder = ['/home/data/subjects/' subject '/raw/'];
        preprocfolder = ['/home/data/subjects/' subject '/preprocessed/']; mkdir(preprocfolder)
        
        sessions = dir([rawfolder '/']); 
        folderindex = struct2cell(sessions); folderindex = logical(cell2mat(folderindex(4,:)));
        sessions = sessions(folderindex);
        sessions = sessions(3:end);
        
        all_T1files = cell(0,1);
        firstT1 = true;
        
        for sess = 1:length(sessions)
                
                T1files = dir([rawfolder '/' sessions(sess).name '/o*' T1string '*.nii.gz']);
                for file = 1:length(T1files)
                    filename = [rawfolder '/' sessions(sess).name '/' T1files(file).name];
                    [~,sizestr] = system(['fslsize ' filename ' -s']);
                    tokens = tokenize(sizestr,' ');
                    datasize = [str2num(tokens{3}) str2num(tokens{5}) str2num(tokens{7}) str2num(tokens{9})];
                    if firstT1
                        T1size = datasize;
                        all_T1files{end+1} = filename;
                        firstT1 = false;
                    else
                        if all(T1size==datasize)
                            all_T1files{end+1} = filename;
                        end
                    end
                end
        end
        
        if length(all_T1files) > 1
            fslviewstring = 'fslview';
            for filenum = 1:length(all_T1files)
                fslviewstring = [fslviewstring ' ' all_T1files{filenum}];
            end
            [~,~] = system(fslviewstring);
            bestnum = input(['Subject ' subject ': input index of best image (counting from bottom of fslview list): ']);
        else
            bestnum = 1;
        end
        copyfile(all_T1files{bestnum},[preprocfolder '/T1.nii.gz']);
    end
end



warning off


if any([copy_rawBOLDdata T1_preprocess T1_segment T2_preprocess surface_registration BOLD_preprocess])

disp('PREPROCESSING')

%% Set up parallel processing
maxworkers = 11;
nworkers = min([subcount maxworkers]);

%processingpool = parpool(nworkers);

%parfor subnum = 1:subcount
for subnum = 1:subcount
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    warning off
    
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
        
        subQCindices = find(strcmp(QCdata(:,1),subject));
        subQC_sessions = cell2mat(QCdata(subQCindices,2));
        
        if copy_rawBOLDdata
            disp(['Subject ' subject ': finding and copying BOLD data']);
            dlmwrite(logfile,'copying raw BOLD data to preprocessing folder','-append','delimiter','')
        end
            
            
        for sess = 1:length(sessions)
            
            dashloc = strfind(sessions(sess).name,'-');
            sessnum = str2num(sessions(sess).name(dashloc+1:end));
            sessQCindices = subQCindices(subQC_sessions==sessnum);
            sessQC_decision = cell2mat(QCdata(sessQCindices,7)); sessQC_decision(isnan(sessQC_decision)) = 1;
            sessQC_acquisitionnum = cell2mat(QCdata(sessQCindices,8)); sessQC_acquisitionnum(isnan(sessQC_acquisitionnum)) = 0;
            
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
                    
                    %check if data has been flagged by manual QC check
                    flaggedasbad = (sessQC_decision(sessQC_acquisitionnum==acquisitionnum)==0);
                    if ~any(sessQC_acquisitionnum==acquisitionnum)
                        flaggedasbad = false;
                    end
                    
                    if correctsize && data_islong && (~flaggedasbad)
                        
                        %keep track of session names, session number for each sequence, and overall session number
                        seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
                        BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
                        sesscounter(end+1) = sess;
                        
                        if copy_rawBOLDdata
                            %load data, trim first few frames
                            data = load_untouch_nii(filename);
                            data.img = data.img(:,:,:,skipframes+1:end);
                            data.hdr.dime.dim(5) = size(data.img,4);
                            
                            %save to preprocessing folder
                            save_untouch_nii(data,[BOLDrunnames{end} '.nii.gz'])
                        end
                        
                    end
                    
                end
            end
            
        end
            
        
        
        
        
        
        %% T1 Preprocessing
        
        
        T1file = 'T1';
        
        if T1_preprocess
            
            cd(preprocfolder)
            
            dlmwrite(logfile,'preprocessing T1','-append','delimiter','')
            
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting, and skull-stripping T1, and registering to MNI']);
            
            %reorient, crop, and bias-correct T1 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T1file ' -o struct --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            copyfile([preprocfolder '/struct.anat/T1_biascorr.nii.gz'],[preprocfolder '/' T1file '_biascorr.nii.gz']);
             T1file = [T1file '_biascorr'];
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['chmod 777 -R ' preprocfolder '/struct.anat']);
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
            ign = rmdir([freesurferfolder '/' subject ],'s');
            cd(preprocfolder)
            
            if logical(systemresult{end,1})
               error('Problem detected with Freesurfer')
            end

        end
        
        
        %% Surface-based registration
        
        if surface_registration
            
            disp(['Subject ' subject ': surface-based registration']);
            dlmwrite(logfile,'conducting surface registration','-append','delimiter','')
            
            %get segmentation masks into nifti format
            systemresult = make_fs_masks(freesurferfolder,normalized_voxdim,[preprocfolder '/T1_2MNI.mat'],template,systemresult);
            copyfile([freesurferfolder '/nusmask/aparc+aseg_cerebralwm.nii.gz'],preprocfolder)
            
            mkdir(fsLRfolder)
            
            prevsize = size(systemresult,1);
            systemresult = PostFreeSurferPipeline_fsavg2fslr_long_singlesub(subject,T1file,freesurferfolder,fsLRfolder,preprocfolder,systemresult);
            if any(cell2mat(systemresult(prevsize+1:end,1)))
                error('Problem detected with surface registration')
            end
            
            
            
            surfacetousefolder = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/'];
            distfolder = [surfacetousefolder '/distances/'];
            mkdir(distfolder)
            
            surfcoordsfileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR.coord.gii'];
            surftopofileL = [surfacetousefolder '/' subject '.L.32k_fs_LR.topo.gii'];
            
            delete([distfolder '/Surface_distances_L_noHEAD.func.gii'])
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -surface-geodesic ' surfcoordsfileL ' ' surftopofileL ' ' distfolder '/Surface_distances_L.func.gii true']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -file-convert -format-convert ASCII ' distfolder '/Surface_distances_L.func.gii']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['awk ''NF > 25'' ' distfolder '/Surface_distances_L.func.gii > ' distfolder '/Surface_distances_L_noHEAD.func.gii']);
            delete([distfolder '/Surface_distances_L.func.gii'])
            
            
            surfcoordsfileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR.coord.gii'];
            surftopofileR = [surfacetousefolder '/' subject '.R.32k_fs_LR.topo.gii'];
            
            delete([distfolder '/Surface_distances_R_noHEAD.func.gii'])
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -surface-geodesic ' surfcoordsfileR ' ' surftopofileR ' ' distfolder '/Surface_distances_R.func.gii true']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['caret_command -file-convert -format-convert ASCII ' distfolder '/Surface_distances_R.func.gii']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['awk ''NF > 25'' ' distfolder '/Surface_distances_R.func.gii > ' distfolder '/Surface_distances_R_noHEAD.func.gii']);
            delete([distfolder '/Surface_distances_R.func.gii'])
            
            
            surffileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR.surf.gii'];
            surfaceareafileL = [surfacetousefolder '/' subject '.L.midthickness.32k_fs_LR_surfaceareas.func.gii'];
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['wb_command -surface-vertex-areas ' surffileL ' ' surfaceareafileL]);
            
            surffileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR.surf.gii'];
            surfaceareafileR = [surfacetousefolder '/' subject '.R.midthickness.32k_fs_LR_surfaceareas.func.gii'];
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['wb_command -surface-vertex-areas ' surffileR ' ' surfaceareafileR]);
            
        end
        
        %% T2 Preprocessing
        
        T2file = 'T2';
        
        if T2_preprocess
            
            cd(preprocfolder)
            
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
        
        
        %% Functional preprocessing
        
        
        
        if BOLD_preprocess
            
            cd(preprocfolder)
            
            dlmwrite(logfile,'preprocessing BOLD data','-append','delimiter','')
            
            %make a file that tracks session number
            sess_tracker = [preprocfolder '/BOLDruns_sessions.txt'];
            delete(sess_tracker);
            fid = fopen(sess_tracker,'at'); %open the output file for writing
            fclose(fid);
            
            brainmask = load_untouch_nii([freesurferfolder '/nusmask/aparc+aseg_brainmask_mask_' num2str(normalized_voxdim) num2str(normalized_voxdim) num2str(normalized_voxdim) '.nii.gz']);
            
            avg_mergestr_native = 'fslmerge -t allBOLDavgs ';
            avg_mergestr_MNI = 'fslmerge -t allBOLDavgs_MNI ';
            
            disp(['Subject ' subject ': slice time correction, motion correction, atlas registration, and intensity normalization']);
            
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
                
                %motion-correct raw BOLD images to first image to get "true" motion plots
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' BOLDfile ' -refvol 0 -plots']);
                
                %delete motion corrected image--we will include this in a single-step atlas transformation below
                delete([BOLDfile '_mcf.nii.gz'])
                
                %slice-time correct BOLD images
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['slicetimer -i ' BOLDfile ' -o ' BOLDfile '_st -r ' num2str(TR) ' --ocustom=' sliceorderfile]);
                
                ign = rmdir([BOLDfile '_st_mcf.mat*'],'s');
                %motion-correct BOLD images to first image to get motion plots
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' BOLDfile '_st -refvol 0 -mats']);
                
                %split slice-time corrected timeseries into separate volumes
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslsplit ' BOLDfile '_st ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf -t']);
                
                %get avg BOLD for this run
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' BOLDfile '_st_mcf -Tmean ' BOLDfile '_st_mcf_avg']);
                
                %Add to string to concatenate avg BOLD images
                avg_mergestr_native = [avg_mergestr_native ' ' BOLDfile '_st_mcf_avg'];
                
                %delete motion corrected image--we will include this in a single-step atlas transformation below
                delete([BOLDfile '_st_mcf.nii.gz'])
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' BOLDfile '_st_mcf_avg ' BOLDfile '_st_mcf_avg_bet -R']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg_bet -ref '  T1file '_bet -dof 12 -omat ' BOLDfile '_avg_2T1_init.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg_bet -ref '  T1file '_bet -dof 12 -cost bbr -wmseg aparc+aseg_cerebralwm -init ' BOLDfile '_avg_2T1_init.mat -omat ' BOLDfile '_avg_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_avg_2MNI.mat -concat T1_2MNI.mat ' BOLDfile '_avg_2T1.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg -interp spline -ref ' template ' -applyisoxfm ' num2str(normalized_voxdim) ' -init ' BOLDfile '_avg_2MNI.mat -out ' BOLDfile '_st_mcf_avg_MNI']);
                
                avg_mergestr_MNI = [avg_mergestr_MNI ' ' BOLDfile '_st_mcf_avg_MNI'];
                
                %get volume numbers
                vols = dir([BOLDfile '_st_mcf.mat/MAT*']);
                fslmergestr = ['fslmerge -t ' BOLDfile '_st_mcf_MNI'];
                %for each volume
                for i = 1:length(vols)
                    %concatenate motion correction transform with run avg to MNI transform
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat -concat  ' BOLDfile '_avg_2MNI.mat ' BOLDfile '_st_mcf.mat/' vols(i).name]);
                    
                    %apply concatenated transform
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf' sprintf('%04i',(i-1)) ' -interp spline -ref ' template ' -applyisoxfm ' num2str(normalized_voxdim) ' -init '  BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat -out ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf_MNI' sprintf('%04i',(i-1))]);
                    
                    %prepare to merge all volumes back into a 4D image
                    fslmergestr = [fslmergestr ' ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf_MNI' sprintf('%04i',(i-1))];
                end
                %merge volumes
                [systemresult{end+1,1},systemresult{end+1,2}] = system(fslmergestr);
                
                 %remove unneeded files
                delete([BOLDfile '_st_mcf.mat/*'])
                rmdir([BOLDfile '_st_mcf.mat'],'s')
                delete([BOLDfile '_st.nii.gz'])
                delete([BOLDfile '_avg_2T1_init.mat'])
                delete([BOLDfile '_st_mcf_avg_bet.nii.gz'])
                
                
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
                dlmwrite(sess_tracker,[BOLDfile '_st_mcf_MNI.nii.gz ' num2str(sesscounter(runnum)) ' ' BOLDfile '_mcf.par'],'-append','delimiter','');%write the data to the output file
                
             end
            
            %merge averages
            [systemresult{end+1,1},systemresult{end+1,2}] = system(avg_mergestr_native);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(avg_mergestr_MNI);

            
            brainmask = [];
            
        end
        
        disp(['Subject ' subject ' COMPLETED preprocessing without error!'])
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


try delete(processingpool); catch; end


end










if any([BOLD_fcprocess cifti_creation cifti_correlation])

    enoughtime = true(length(subjects),length(functional_sequences));
    
    disp('POSTPROCESSING')
    
for subnum = 1:subcount
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    warning off
    
    try
        
        %% Set up folders
        
        
        rawfolder = ['/home/data/subjects/' subject '/raw/'];
        preprocfolder = ['/home/data/subjects/' subject '/preprocessed/'];
        freesurferfolder = ['/home/data/subjects/' subject '/freesurfer/']; 
        fsLRfolder = ['/home/data/subjects/' subject '/fs_LR/']; 
        fc_processfolder = ['/home/data/subjects/' subject '/fc_processed/']; 
        ciftifolder = ['/home/data/subjects/' subject '/cifti/']; 
        correlationfolder = ['/home/data/subjects/' subject '/connectome/'];
        infomapfolder = ['/home/data/subjects/' subject '/infomap/'];
        logfolder = ['/home/data/subjects/' subject '/processing_logs/']; 
        
        logfile2 = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile2,' '); logfile2 = [tokens{1} '_' tokens{2}];
        
        %% Fc-processing
        
        
        if BOLD_fcprocess
            
            mkdir(fc_processfolder)
            
            dlmwrite(logfile2,'fc-processing BOLD data','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                if fcprocess_sequences{seq}
                    
                    disp(['Subject ' subject ', ' functional_sequences{seq} ': fc-processing']);
                    FC_Process(subject,['/home/data/subjects/' subject '/preprocessed/BOLDruns_sessions.txt'],FDthresh,functional_sequence_TRs{seq},['/home/data/subjects/' subject '/fc_processed/'],functional_sequences{seq})
                    
                end
            end
            
        end
        
        
        for seq = 1:length(functional_sequences)
            tmask = load([fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt']);
            if (nnz(tmask) .* functional_sequence_TRs{seq} ./ 60) < mintime
                enoughtime(subnum,seq) = false;
                disp(['Subject ' subject ' does not have at least ' num2str(mintime) ' minutes of ' functional_sequences{seq} ' data! Processing will not continue.'])
            end
        end
        
        
        
        %% Cifti creation
        
        if cifti_creation
            
            mkdir(ciftifolder)
            
            dlmwrite(logfile2,'creating ciftis','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                if enoughtime(subnum,seq)
                    
                    disp(['Subject ' subject ', ' functional_sequences{seq} ': mapping functional data to cortical surface']);
                    
                    if fcprocess_sequences{seq}
                        
                        systemresult = post_fc_processing_batch_singlesub(subject,[fc_processfolder '/' functional_sequences{seq} '_fc_processed_tmasked.nii.gz'],ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt'],[preprocfolder '/BOLDruns_sessions.txt'],[fc_processfolder '/' functional_sequences{seq} '_runs_sessions.txt'],'/home/data/scripts/Cifti_creation/post_fc_processing_batch_params_singlesub.m',systemresult);
                        
                    else
                        
                        BOLDmergedfile = [preprocfolder '/' functional_sequences{seq} '_merged'];
                        
                        fslmergestr = ['fslmerge -t ' BOLDmergedfile ' '];
                        
                        [BOLDrunnames,~,~] = textread([preprocfolder '/BOLDruns_sessions.txt'],'%s%s%s');
                        
                        for runnum = 1:length(BOLDrunnames)
                            [~,BOLDfilename,~] = fileparts(BOLDrunnames{runnum});
                            if strcmp(BOLDfilename(1:length(functional_sequences{seq})),functional_sequences{seq});
                                fslmergestr = [fslmergestr BOLDrunnames{runnum} '_st_mcf_MNI.nii.gz '];
                            end
                        end
                        [systemresult{end+1,1},systemresult{end+1,2}] = system(fslmergestr);
                        
                        systemresult = post_fc_processing_batch_singlesub(subject,BOLDmergedfile,ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[],[],[],'/home/data/scripts/Cifti_creation/post_fc_processing_batch_params_singlesub.m',systemresult);
                        
                        delete(BOLDmergedfile)
                        
                    end
                    if systemresult{end,1}
                        error('Problem detected with Cifti creation')
                    end
                end
            end
            
            
            
        end
        
        
        
        
        %% Cifti correlation
        
        if cifti_correlation
           
            mkdir(correlationfolder)
            dlmwrite(logfile2,'correlating ciftis and conducting geodesic distance calculation','-append','delimiter','')
            
            
            for seq = 1:length(functional_sequences)
                
                if enoughtime(subnum,seq)
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': correlating all timecourses']);
                 
                 ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
                 
                data = ft_read_cifti_mod(ciftifile);
                
                corr = paircorr_mod(data.data');
                corr = FisherTransform(corr);
                corr(isnan(corr)) = 0;
                
                data.dimord = 'pos_pos';
                data.data = corr;
                clear corr
                
                ft_write_cifti_mod([correlationfolder '/' functional_sequences{seq} '_corr.dconn.nii'],data)
                
                data = [];
                
                end
                
            end
            
            if any(enoughtime(subnum,:))
            
            disp(['Subject ' subject ': getting point-to-point geodesic and euclidean distances'])
            
            surfacetousefolder = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/'];
            distfolder = [surfacetousefolder '/distances/'];
            mkdir(distfolder)
            
            distancesLfile = [distfolder '/Surface_distances_L.mat'];
            if ~exist(distancesLfile)
            geo_distances = load([distfolder '/Surface_distances_L_noHEAD.func.gii']); 
            geo_distances = single(geo_distances); geo_distances(:,1) = [];
            save(distancesLfile,'geo_distances','-v7.3');%parsave(distancesLfile,geo_distances,'-v7.3')
            clear geo_distances
            delete([distfolder '/Surface_distances_L_noHEAD.func.gii'])
            end
            
            distancesRfile = [distfolder '/Surface_distances_R.mat'];
            if ~exist(distancesRfile)
            geo_distances = load([distfolder '/Surface_distances_R_noHEAD.func.gii']);
            geo_distances = single(geo_distances); geo_distances(:,1) = [];
            save(distancesRfile,'geo_distances','-v7.3');%parsave(distancesRfile,geo_distances,'-v7.3')
            clear geo_distances
            delete([distfolder '/Surface_distances_R_noHEAD.func.gii'])
            end
            
            mkdir([ciftifolder '/distances/'])
            Make_distmat_32k_fsLR_singlesub(subject,ciftifile,[ciftifolder '/distances/normalwall_distmat.mat'])
            
            end
            
            
        end
        

                
                
                
            
            
        %%
        
        
        disp(['Subject ' subject ' COMPLETED post-processing without error!'])
        dlmwrite(logfile2,'all processing complete.','-append','delimiter','')
        
    catch errorinfo
        
        dlmwrite(logfile2,'ERROR:','-append','delimiter','')
        
        errormessage = [errorinfo.stack(1).file ', line ' num2str(errorinfo.stack(1).line) ': ' errorinfo.message];
        
        probleminds = find(cell2mat(systemresult(:,1)));
        for i = 1:length(probleminds)
            dlmwrite(logfile2,systemresult{probleminds(i),2},'-append','delimiter','');
        end
        dlmwrite(logfile2,errormessage,'-append','delimiter','');
        
        disp(['Subject ' subject ' FAILED to complete post-processing! Check error log ' logfile2])
        
    end
    
    
end

end

