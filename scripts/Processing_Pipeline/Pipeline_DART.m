%% Parameters

template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz';

T1string = 'prage';
T2string = 'T2';

normalized_voxdim = 3;

functional_sequences = {'RSFC'};%,'Study','Restudy'};%{'restingstate'};

functional_sequence_TRs = {2.5,2.5,2.5};

functional_sequence_rawsize = {[64 64 32],[64 64 32],[64 64 32]};

fcprocess_sequences = {true,false,false};

skipframes = 5;

FDthresh = .2;

mintime = 0; % minutes

infomapthresholds = [.005 : .0025 : .05];

sliceorderfile = '/home/data/scripts/Resources/sliceorder_interleaved_even.txt';

%smallwall_name = 'DART';




%% Steps to run

choose_T1 = false;
copy_rawdata = false;
T1_preprocess = false;
T1_segment = false;
T2_preprocess = false;
surface_registration = false;
BOLD_preprocess = false;
%smallwall_creation = true;
BOLD_fcprocess = false;
cifti_creation = true;
cifti_correlation = false;
infomap = false;


%% Get subjects

%subjects = '/home/data/subjects/DART.txt';
%subjects = '/home/data/subjects/processing_list_allcomplete.txt';
%subjects = '/home/data/subjects/list_forJDP.txt';
subjects = {'vc35498'};

if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end

subcount = length(subjects);

seq_counter = cell(subcount,1);

warning off

%% Choose best T1 image

if choose_T1
    disp('Copying T1 images')
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


%% Set up parallel processing
warning off

maxworkers = 20;

if any([copy_rawdata T1_preprocess T1_segment T2_preprocess surface_registration BOLD_preprocess])

disp('PREPROCESSING')

nworkers = min([subcount maxworkers]);
    
processingpool = parpool(nworkers);

parfor subnum = 1:subcount
%for subnum = 1:subcount
    
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
        
        if copy_rawdata
            disp(['Subject ' subject ': finding and copying data']);
            dlmwrite(logfile,'copying raw BOLD data to preprocessing folder','-append','delimiter','')
            
            foundfieldmaps = false;
            
            for sess = 1:length(sessions)
                
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
               
                if ~foundfieldmaps
                
                fieldmapstr = 'fieldmap';
                fieldmaps = dir([rawfolder '/' sessions(sess).name '/*' fieldmapstr '*.nii.gz']);
                if length(fieldmaps)==2
                    foundfieldmaps = true;
                    copyfile([rawfolder '/' sessions(sess).name '/' fieldmaps(1).name],[preprocfolder '/fieldmap_mag.nii.gz'])
                    copyfile([rawfolder '/' sessions(sess).name '/' fieldmaps(2).name],[preprocfolder '/fieldmap_phase.nii.gz'])
                end
                end
                
                
            end
            
            
            
            
            
            
            
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
        
        
        T1file = 'T1';
        
        if T1_preprocess
            
            cd(preprocfolder)
            
            dlmwrite(logfile,'preprocessing T1','-append','delimiter','')
            
            disp(['Subject ' subject ': reorienting, cropping, bias-correcting, and skull-stripping T1, and registering to MNI']);
            
            %reorient, crop, and bias-correct T1 image
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['fsl_anat -i ' T1file ' -o struct --clobber --noreg --nononlinreg --noseg --nosubcortseg']);
            copyfile([preprocfolder '/struct.anat/T1_biascorr.nii.gz'],[preprocfolder '/' T1file '_biascorr.nii.gz'])
            T1file = [T1file '_biascorr'];
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['chmod 777 -R ' preprocfolder '/struct.anat']);
            delete([preprocfolder '/struct.anat/*']);
            [~] = rmdir([preprocfolder '/struct.anat'],'s');
            
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
        
        
        
        T2file = 'T2_avg';
        
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
                
            
            %prepare fieldmap
            [systemresult{end+1,1},systemresult{end+1,2}] = system('bet fieldmap_mag fieldmap_mag_bet -R');
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in fieldmap_mag_bet -ref ' T1file '_bet -dof 12 -omat fieldmap_2T1.mat']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system('fslmaths fieldmap_mag_bet -ero fieldmap_mag_bet');
            [systemresult{end+1,1},systemresult{end+1,2}] = system('fslmaths fieldmap_mag_bet -ero fieldmap_mag_bet');
            [systemresult{end+1,1},systemresult{end+1,2}] = system('fslmaths fieldmap_mag_bet -bin fieldmap_mag_bet_mask');
            [systemresult{end+1,1},systemresult{end+1,2}] = system('fsl_prepare_fieldmap SIEMENS fieldmap_phase fieldmap_mag_bet fieldmap_rads 2.46');
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in fieldmap_rads -ref ' T1file '_bet -applyxfm -init fieldmap_2T1.mat -out fieldmap_rads_T1space']);
            [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in fieldmap_mag_bet_mask -ref ' T1file '_bet -applyxfm -init fieldmap_2T1.mat -out fieldmap_mag_bet_mask_T1space']);
                
                
            for runnum = 1:length(BOLDrunnames)
                BOLDfile = BOLDrunnames{runnum};
                
                %get short filename
                [~,BOLDfilename,~] = fileparts(BOLDfile);
                
                %get TR for this run
                for i = 1:length(functional_sequences)
                    if (length(BOLDfilename) >= length(functional_sequences{i})) && (strcmp(BOLDfilename(1:length(functional_sequences{i})),functional_sequences{i}));
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
                

                %atlas registration
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['bet ' BOLDfile '_st_mcf_avg ' BOLDfile '_st_mcf_avg_bet -R']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg_bet -ref '  T1file '_bet -dof 12 -omat ' BOLDfile '_avg_2T1_init.mat']);
                [systemresult{end+1,1},systemresult{end+1,2}] = system(['flirt -in ' BOLDfile '_st_mcf_avg_bet -ref '  T1file '_bet -dof 12 -cost bbr -wmseg aparc+aseg_cerebralwm -init ' BOLDfile '_avg_2T1_init.mat -omat ' BOLDfile '_avg_2T1.mat -schedule ${FSLDIR}/etc/flirtsch/bbr.sch -fieldmap fieldmap_rads_T1space -fieldmapmask fieldmap_mag_bet_mask_T1space -pedir -2 -echospacing .00059']);
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
        
        errormessage = [errorinfo.stack(1).file ', line ' num2str(errorinfo.stack(1).line) ': ' errorinfo.message];
        
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




% if smallwall_creation
%     for subnum = 1:length(subjects)
%         thissub_L_white = gifti(['/home/data/subjects/' subjects{subnum} '/fs_LR/MNI/fsaverage_LR32k/' subjects{subnum} '.L.white.32k_fs_LR.surf.gii']);
%         thissub_L_pial = gifti(['/home/data/subjects/' subjects{subnum} '/fs_LR/MNI/fsaverage_LR32k/' subjects{subnum} '.L.pial.32k_fs_LR.surf.gii']);
%         sub_walls_L_all(:,subnum) = all(thissub_L_white.vertices==thissub_L_pial.vertices,2);
%         
%         thissub_R_white = gifti(['/home/data/subjects/' subjects{subnum} '/fs_LR/MNI/fsaverage_LR32k/' subjects{subnum} '.R.white.32k_fs_LR.surf.gii']);
%         thissub_R_pial = gifti(['/home/data/subjects/' subjects{subnum} '/fs_LR/MNI/fsaverage_LR32k/' subjects{subnum} '.R.pial.32k_fs_LR.surf.gii']);
%         sub_walls_R_all(:,subnum) = all(thissub_R_white.vertices==thissub_R_pial.vertices,2);
%     end
%     walls_L = any(sub_walls_L_all,2);
%     walls_R = any(sub_walls_R_all,2);
%     save(gifti(single(walls_L)),['/home/data/scripts/Resources/cifti_masks/' smallwall_name 'L.commonverts.32k_fs_LR.shape.gii'])
%     save(gifti(single(walls_R)),['/home/data/scripts/Resources/cifti_masks/' smallwall_name 'R.commonverts.32k_fs_LR.shape.gii'])
% end





if any([BOLD_fcprocess cifti_creation cifti_correlation infomap])

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
            if fcprocess_sequences{seq}
                tmask = load([fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt']);
                if (nnz(tmask) .* functional_sequence_TRs{seq} / 60) < mintime
                    enoughtime(subnum,seq) = false;
                    disp(['Subject ' subject ' does not have at least ' num2str(mintime) ' minutes of ' functional_sequences{seq} ' data! Processing will not continue.'])
                end
            end
        end
        
        
        
        %% Cifti creation
        
        if cifti_creation
            
            mkdir(ciftifolder)
            
            dlmwrite(logfile2,'creating ciftis','-append','delimiter','')
            
            
            
            medial_masks = {'/home/data/scripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii','/home/data/scripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';...
                '/home/data/scripts/Resources/cifti_masks/L.atlasroi_erode3.32k_fs_LR.shape.gii','/home/data/scripts/Resources/cifti_masks/R.atlasroi_erode3.32k_fs_LR.shape.gii'};
            
            for seq = 1:length(functional_sequences)
                
                if enoughtime(subnum,seq)
                    
                    disp(['Subject ' subject ', ' functional_sequences{seq} ': mapping functional data to cortical surface']);
                    
                    if fcprocess_sequences{seq}
                        
                        %systemresult = post_fc_processing_batch_singlesub(subject,[fc_processfolder '/' functional_sequences{seq} '_fc_processed_tmasked.nii.gz'],ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt'],[preprocfolder '/BOLDruns_sessions.txt'],[fc_processfolder '/' functional_sequences{seq} '_runs_sessions.txt'],medial_masks,systemresult);
                        systemresult = post_fc_processing_batch_singlesub_nosmooth(subject,[fc_processfolder '/' functional_sequences{seq} '_fc_processed_tmasked.nii.gz'],ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[fc_processfolder '/' functional_sequences{seq} '_all_tmask.txt'],[preprocfolder '/BOLDruns_sessions.txt'],[fc_processfolder '/' functional_sequences{seq} '_runs_sessions.txt'],medial_masks,systemresult);
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
                        
                        systemresult = post_fc_processing_batch_singlesub(subject,BOLDmergedfile,ciftifolder,[fsLRfolder '/MNI/'],functional_sequence_TRs{seq},functional_sequences{seq},[],[],[],medial_masks,systemresult);
                        
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
        
        
        
        %% Infomap
        
        if infomap
            
            minsize = 400;
            xdist = 30;
            
            mkdir(infomapfolder)
            
            dlmwrite(logfile2,'running infomap','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                if enoughtime(subnum,seq)
                
                outfolder = [infomapfolder '/' functional_sequences{seq} '/'];
                mkdir(outfolder)
            
                disp(['Subject ' subject ', ' functional_sequences{seq} ': running infomap']);
                
                corrfile = [correlationfolder '/' functional_sequences{seq} '_corr.dconn.nii'];
                
                dmatfile = [ciftifolder '/distances/normalwall_distmat_uint8.mat'];
                
                Run_Infomap(corrfile, dmatfile, xdist, infomapthresholds, 0, outfolder);
                
                
                simplified=modify_clrfile('simplify',[outfolder '/rawassn.txt'],minsize);
                save([outfolder '/rawassn_minsize' num2str(minsize) '.txt'],simplified);
                regularized =  regularize(simplified,[outfolder '/rawassn_minsize' num2str(minsize) '_regularized.txt']);
                
                outcifti = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
                outcifti.data = regularized;
                outcifti.dimord = 'scalar_pos';
                for i = 1:length(infomapthresholds)
                    outcifti.mapname{1,i} = ['minsize=' num2str(minsize) '; xdist=' num2str(xdist) ';thresh=' num2str(infomapthresholds(i))];
                end
                ft_write_cifti_mod([outfolder '/rawassn_minsize' num2str(minsize) '_regularized.dscalar.nii'],outcifti);
                
                end
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

