%% Parameters

template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz';

normalized_voxdim = 3;

functional_sequences = {'resting_state'};

functional_sequence_TRs = {2.2};

fcprocess_sequences = {false};

skipframes = 5;




%% Steps to run

convert_dicom_data = false;
copy_rawdata = true;
T1_preprocess = true;
T1_segment = false;
T2_preprocess = false;
surface_registration = false;
BOLD_preprocess = true;
BOLD_fcprocess = true;
cifti_creation = true;


%% Get subjects

subjects = {'SMN'};

if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end


%% Set up parallel processing

seq_counter = cell(nsubs,1);

maxworkers = 24;

nsubs = length(subjects);
nworkers = min([nsubs maxworkers]);

processingpool = parpool(nworkers);

parfor subnum = 1:nsubs
%for subnum = 1:nsubs
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    warning off
    errorfile = ['/data/subjects/' subject '/processing_errors.txt'];
    delete(errorfile);
    
    try
        
        %% Set up folders
        
        
        rawfolder = ['/data/subjects/' subject '/raw/'];
        preprocfolder = ['/data/subjects/' subject '/preprocessed/']; mkdir(preprocfolder)
        freesurferfolder = ['/data/subjects/' subject '/freesurfer/']; mkdir(freesurferfolder)
        fsLRfolder = ['/data/subjects/' subject '/fs_LR/']; mkdir(fsLRfolder)
        fc_processfolder = ['/data/subjects/' subject '/fc_processed/']; mkdir(fc_processfolder)
        ciftifolder = ['/data/subjects/' subject '/cifti/']; mkdir(ciftifolder)
        
        %% Convert dicom data to nifti and anonymize
        
        if convert_dicom_data
           
           disp(['Subject ' subject ': converting and anonymizing data']);
           
           [systemresult{end+1,1},systemresult{end+1,2}] = system(['python /data/scripts/Processing_Pipeline/anon_dcn_files -b ' rawfolder ' -o ' rawfolder]);
           
            
        end
        
        %% Find nifti data in raw folder and move to preprocessed folder
        
        sesscounter = [];
        
        sessions = dir([rawfolder '/*_anon']);
        
        for seq = 1:length(functional_sequences)
            seq_counter{subnum}{seq} = 0;
        end
        BOLDrunnames = cell(0,1);
        
        if copy_rawdata
            
            T1_concatstr = ['fslmerge -t ' preprocfolder '/T1_concat'];
            T2_concatstr = ['fslmerge -t ' preprocfolder '/T2_concat'];
            
            for sess = 1:length(sessions)
                
                %Find all T1s and prepare to align and average them
                T1folders = dir([rawfolder '/' sessions(sess).name '/NIFTI/*T1*']);
                for folder = 1:length(T1folders);
                    filename = dir([rawfolder '/' sessions(sess).name '/NIFTI/' T1folders(folder).name '/*.nii.gz']);
                    T1_concatstr = [T1_concatstr ' ' rawfolder '/' sessions(sess).name '/NIFTI/' T1folders(folder).name '/' filename(1).name];
                end
                
                %Find all T2s and prepare to align and average them
                T2folders = dir([rawfolder '/' sessions(sess).name '/NIFTI/*T2*']);
                for folder = 1:length(T2folders);
                    filename = dir([rawfolder '/' sessions(sess).name '/NIFTI/' T2folders(folder).name '/*.nii.gz']);
                    T2_concatstr = [T2_concatstr ' ' rawfolder '/' sessions(sess).name '/NIFTI/' T2folders(folder).name '/' filename(1).name];
                end
                
                %Find all BOLD runs specified
                disp(['Subject ' subject ': finding and copying BOLD runs']);
                                                                
                for seq = 1:length(functional_sequences)
                    thisseq_folders = dir([rawfolder '/' sessions(sess).name '/NIFTI/*' functional_sequences{seq} '*']);
                    for folder = 1:length(thisseq_folders);
                        filename = dir([rawfolder '/' sessions(sess).name '/NIFTI/' thisseq_folders(folder).name '/*.nii.gz']);
                        
                        %keep track of session names, session number for each sequence, and overall session number
                        seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
                        BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
                        sesscounter(end+1) = sess;
                        
                        %load data, trim first few frames, and save to preprocessing folder
                        data = load_untouch_nii([rawfolder '/' sessions(sess).name '/NIFTI/' thisseq_folders(folder).name '/' filename(1).name]);
                        data.img = data.img(:,:,:,skipframes+1:end);
                        data.dime.dim(5) = size(data.img,4);
                        save_untouch_nii(data,[BOLDrunnames{end} '.nii.gz'])
                        
                    end
                end
                
            end
            
            %T1s: concatenate, align, and average
            disp(['Subject ' subject ': aligning and averaging T1s']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(T1_concatstr);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' preprocfolder '/T1_concat -refvol 0']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' preprocfolder '/T1_concat_mcf -Tmean ' preprocfolder 'T1_avg']);
            delete([preprocfolder '/T1_concat.nii.gz'])
            delete([preprocfolder '/T1_concat_mcf.nii.gz'])
            
            %T2s: concatenate, align, and average
            disp(['Subject ' subject ': aligning and averaging T2s']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(T2_concatstr);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' preprocfolder '/T2_concat -refvol 0']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' preprocfolder '/T2_concat_mcf -Tmean ' preprocfolder 'T2_avg']);
            delete([preprocfolder '/T2_concat.nii.gz'])
            delete([preprocfolder '/T2_concat_mcf.nii.gz'])
            
        else
            %keep track of functional session names, session number for each sequence, and overall session number
            for sess = 1:length(sessions)
                for seq = 1:length(functional_sequences)
                    thisseq_folders = dir([rawfolder '/' sessions(sess).name '/NIFTI/*' functional_sequences{seq} '*']);
                    for folder = 1:length(thisseq_folders);
                        seq_counter{subnum}{seq} = seq_counter{subnum}{seq} +1;
                        BOLDrunnames(end+1) = {[preprocfolder '/' functional_sequences{seq} '_' num2str(seq_counter{subnum}{seq})]};
                        sesscounter(end+1) = sess;
                    end
                end
            end
            
            
        end
        
        
        
        
        %% T1 Preprocessing
        cd(preprocfolder)
        
        T1file = 'T1_avg';
        
        if T1_preprocess
            
            %reorient, crop, and bias-correct T1 image
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting T1']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T1file ' -o struct --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            copyfile([preprocfolder '/struct.anat/T1_biascorr.nii.gz'],[preprocfolder '/' T1file '_biascorr.nii.gz'])
            T1file = [T1file '_biascorr'];
            delete([preprocfolder '/struct.anat/*']);
            rmdir([preprocfolder '/struct.anat']);
            
            %brain-extract T1 image
            disp(['Subject ' subject ': skull-stripping T1']); pause(.1)
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' T1file ' ' T1file '_bet -R -g -.1']);
            
            %register T1 to MNI
            disp(['Subject ' subject ': registering T1 to MNI']); pause(.1)
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T1file '_bet -interp spline -dof 12 -ref ' template ' -omat T1_2MNI.mat -out ' T1file '_MNI']);
            
        else
            
            %track name
            T1file = [T1file '_biascorr'];
            
        end
        
        
        %% T1 Segmentation
        
        if T1_segment
            
            disp(['Subject ' subject ': segmenting T1']);
            
            %copy bias-corrected T1 to freesurfer folder
            copyfile([T1file '.nii.gz'],[freesurferfolder '/' T1file '.nii.gz'])
            cd(freesurferfolder)
            
            %transform T1 to .mgz format
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['recon-all -subjid ' subject ' -sd ' freesurferfolder ' -i ' freesurferfolder '/' subject '_' T1file '.nii.gz']);
            %segment T1
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['recon-all -subjid ' subject ' -sd ' freesurferfolder ' -all']);
            
            %move files around
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['mv ' freesurferfolder '/' subject '/* ' freesurferfolder]);
            rmdir([freesurferfolder '/' subject ])
            cd(preprocfolder)
            
            %get segmentation masks into nifti format
            systemresult = make_fs_masks(freesurferfolder,normalized_voxdim,[preprocfolder '/T1_2MNI.mat'],template,systemresult);
            copyfile([freesurferfolder '/nusmask/aparc+aseg_cerebralwm.nii.gz'],preprocfolder)
            
        end
        
        
        %% T2 Preprocessing
        
        cd(preprocfolder)
        
        T2file = 'T2_avg';
        
        if T2_preprocess
            
            %reorient, crop, and bias-correct T2 image
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting T2']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T2file ' -o struct -t T2 --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            copyfile([preprocfolder '/struct.anat/T2_biascorr.nii.gz'],[preprocfolder '/' T2file '_biascorr.nii.gz'])
            T2file = [T2file '_biascorr'];
            
            %brain-extract T2 image
            %disp(['Subject ' subject ': skull-stripping T2'])
            %[status(end+1) result{end+1}] = system(['bet ' T2file ' ' T2file '_bet -R -g -.1']);
            
            %register T2 to T1, and concatenate with T1 to MNI
            %disp(['Subject ' subject ': registering T2 to T1, and through to MNI'])
            %OPTION 1
            %[status(end+1) result{end+1}] = system(['flirt -in ' T2file ' -ref ' T1file ' -interp spline -dof 6 -omat T2_2T1.mat -out T2_T1space']);
            
            
            %OPTION 2
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -dof 6 -omat T2_2T1_init.mat']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' T2file ' -ref ' T1file ' -interp spline -dof 6 -cost bbr -wmseg aparc+aseg_cerebralwm -init T2_2T1_init.mat -omat T2_2T1.mat -out T2_T1space -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
            
            [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat T1_2T2.mat -inverse T2_2T1.mat');
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in aparc+aseg_cerebralwm -ref ' T2file ' -applyxfm -init T1_2T2.mat -interp nearestneighbour -out aparc+aseg_cerebralwm_T2space']);
            
            [systemresult{end+1,1},systemresult{end+1,2}] = system('convert_xfm -omat T2_2MNI.mat -concat T1_2MNI.mat T2_2T1.mat');
            
        else
            
            %track T2 name
            T2file = [T2file '_biascorr'];
            
        end
        
        
        %% Surface-based registration
        
        if surface_registration
            
            disp(['Subject ' subject ': surface-based registration']);
            systemresult = PostFreeSurferPipeline_fsavg2fslr_long_singlesub(subject,T1file,systemresult);
            
        end
        
        
        %% Functional preprocessing
        
        cd(preprocfolder)
        
        if BOLD_preprocess
            
            %make a file that tracks session number
            sess_tracker = [preprocfolder '/BOLDruns_sessions.txt'];
            delete(sess_tracker);
            fid = fopen(sess_tracker,'at'); %open the output file for writing
            fclose(fid);
            
            brainmask = load_untouch_nii([freesurferfolder '/nusmask/aparc+aseg_brainmask_mask_' num2str(normalized_voxdim) num2str(normalized_voxdim) num2str(normalized_voxdim) '.nii.gz']);
            
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
                disp(['Subject ' subject ', ' BOLDfilename ': slice time correction']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['slicetimer -i ' BOLDfile ' -o ' BOLDfile '_st -r ' num2str(TR) ' --ocustom=/data/scripts/Resources/sliceorder.txt']);
                
                %get first (reference) image
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslselectvols -i ' BOLDfile '_st -o ' BOLDfile '_st_vol1 --vols=0']);
                
                %motion-correct BOLD images to first image and save transformation parameters
                disp(['Subject ' subject ', ' BOLDfilename ': motion correction']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['mcflirt -in ' BOLDfile '_st -refvol 0 -mats -plots']);
                
                %delete motion corrected image--we will include this in a single-step atlas transformation below
                delete([BOLDfile '_st_mcf.nii.gz'])
                
                %split slice-time corrected timeseries into separate volumes
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['fslsplit ' BOLDfile '_st ' BOLDfile '_st_mcf.mat/' BOLDfilename '_st_mcf -t']);
                
                
                %register BOLD to T2 and concatenate through to MNI
                disp(['Subject ' subject ', ' BOLDfilename ': register BOLD to T2 and through to MNI']);
                
                %     %OPTION 1
                %     [status(end+1) result{end+1}] = system(['flirt -in ' BOLDfile '_st_mcf -ref '  T2file ' -omat ' BOLDfile '_2T2.mat']);
                %
                %     %OPTION 2
                %     %[status(end+1) result{end+1}] = system(['flirt -in ' BOLDfile '_st_mcf -ref '  T2file '_bet -dof 6 -omat ' BOLDfile '_2T2_init.mat']);
                %     %[status(end+1) result{end+1}] = system(['flirt -in ' BOLDfile '_st_mcf -ref '  T2file '_bet -dof 6 -cost bbr -wmseg aparc+aseg_cerebralwm_T2space -init ' BOLDfile '_2T2_init.mat -omat ' BOLDfile '_2T2.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
                %
                %     [status(end+1) result{end+1}] = system(['convert_xfm -omat ' BOLDfile '_2MNI.mat -concat T2_2MNI.mat ' BOLDfile '_2T2.mat']);
                
                
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_vol1 -ref '  T1file '_bet -dof 6 -omat ' BOLDfile '_2T1_init.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_vol1 -ref '  T1file '_bet -dof 6 -cost bbr -wmseg aparc+aseg_cerebralwm -init ' BOLDfile '_2T1_init.mat -omat ' BOLDfile '_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_2MNI.mat -concat T1_2MNI.mat ' BOLDfile '_2T1.mat']);
                
                disp(['Subject ' subject ', ' BOLDfilename ': apply computed MNI registration']);
                %get volume numbers
                vols = dir([BOLDfile '_st_mcf.mat/MAT*']);
                fslmergestr = ['fslmerge -t ' BOLDfile '_st_mcf_MNI'];
                %for each volume
                for i = 1:length(vols)
                    %concatenate motion correction transform with reference volume to MNI transform
                    [systemresult{end+1,1},systemresult{end+1,2}] = system(['convert_xfm -omat ' BOLDfile '_st_mcf.mat/' vols(i).name '_MNI.mat -concat ' BOLDfile '_st_mcf.mat/' vols(i).name ' ' BOLDfile '_2MNI.mat']);
                    
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
                delete([BOLDfile '_st_vol1.nii.gz'])
                
                
                %Intensity normalization to mode within-brain value of 1000
                disp(['Subject ' subject ', ' BOLDfilename ': intensity normalization to mode 1000']);
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
            brainmask = [];
            
        end
        
        
        %% Fc-processing
        
        
        if BOLD_fcprocess
            
            for seq = 1:length(functional_sequences)
                if fcprocess_sequences{seq}
                    
                    disp(['Subject ' subject ', ' functional_sequences{seq} ': fc-processing']);
                    %FCPROCESS
                    
                end
            end
            
        end
        
        %% Cifti creation
        
        if cifti_creation
            
            for seq = 1:length(functional_sequences)
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': mapping functional data to cortical surface']);
                
                if fcprocess_sequences{seq}
                    
                    systemresult = post_fc_processing_batch_singlesub(subject,voldata_tomap,functional_sequence_TRs{seq},functional_sequences{seq},tmaskfile,preproc_datalist,'/data/scripts/Cifti_creation/post_fc_processing_batch_params_singlesub.m',systemresult);
                    
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
                    
                    systemresult = post_fc_processing_batch_singlesub(subject,BOLDmergedfile,functional_sequence_TRs{seq},functional_sequences{seq},[],[],'/data/scripts/Cifti_creation/post_fc_processing_batch_params_singlesub.m',systemresult);
                    
                    delete(BOLDmergedfile)
                    
                end
            end
            
        end
        
        
        disp(['Subject ' subject ' completed all processing steps without error!'])
        
    catch errorinfo
        
        %Create error tracking file to record issues with processing
        fid = fopen(errorfile,'at');
        fclose(fid);
        
        errormessage = [errorinfo.stack.file ', line ' num2str(errorinfo.stack.line) ': ' errorinfo.message];
        
        probleminds = find(cell2mat(systemresult(:,1)));
        for i = 1:length(probleminds)
            dlmwrite(errorfile,systemresult{probleminds(i),2},'-append','delimiter','');
        end
        dlmwrite(errorfile,errormessage,'-append','delimiter','');
        
        disp(['Subject ' subject ' FAILED to complete processing! Check error log ' errorfile])
        
    end
    
    
end

delete(processingpool)