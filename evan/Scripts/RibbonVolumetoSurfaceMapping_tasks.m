outputdir = '/data/cn4/evan/Task_parcellation/Manytasks_subregister/';


lists = {'/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list'}%,...
%'/data/cn4/evan/Occipitotemporal/Lexicalization/n28_D1_glm_event_cln_gen_anal_wt4.list',...
%     '/data/cn4/evan/Occipitotemporal/Orthography/Adults.glm_list',...
%     '/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt',...
%     '/data/cn4/evan/Task_parcellation/SteveSubjects/Subjectlist_withanat.txt',...'/data/cn3/steven/NP847/glmsLists/28subs_16tps_Prim_list_of_glms.glm_list',...
%     '/data/cn4/maital/reward/REPREPROCESS/33Ss_NotScrubbed/2x3x5_list_of_glms.glm_list',...
%     '/data/cn4/maital/Errors/2TaskGroups/ACRN.glm_list',...
%     '/data/cn3/joe/SustainedTaskLoad/FloatSustainedOnly_Scrubbed_f06_28Subjects.list',...
%     '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...    %'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
%     '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};

directories = {    '/data/cn4/evan/Occipitotemporal/Katie/SingleSubjects/'}%,...
%'/data/cn4/evan/Occipitotemporal/Lexicalization/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/Orthography/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/VisualAttention/SingleSubjects/',...
%     '/data/cn4/evan/Task_parcellation/SteveSubjects/OrigData/',...
%     '/data/cn4/evan/Occipitotemporal/Reward/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/ACRN/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/SustainedTaskLoad/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Trials/SingleSubjects/',...    %'/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/Joe3Tasks_IndConditions/SingleSubjects/'};%,...

extranamestuffstrings = {%'_avgtc_' ...
    %'_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    [],...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...    %'_avgtc_' ...
    '_avgtc_'};

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';

mkdir([outputdir '/goodvoxels'])
        

for listnum = 1:length(lists)
    list = lists{listnum};
    directory = directories{listnum};
    extranamestuff = extranamestuffstrings{listnum};
    if isempty(extranamestuff)
        stevedata = 1;
    else
        stevedata = 0;
    end
    
    if ~isempty(strfind(list,'VisualSwitch'))
        jessdata = 1;
    else
        jessdata = 0;
    end
    
    [front entries] = textread(list,'%s%s','delimiter',':');
    
    %Figure out the number of subjects, the glm files, and the T4 files from
    %the GLM list file
    numsubjects = str2num(entries{2});
    glmfiles = entries(3:2+numsubjects);
    T4files = entries(4+numsubjects:end);
    
    %This part accounts for different GLM list styles
    for i = 1:length(glmfiles)
        if strcmp(glmfiles{i},'')
            glmfiles{i} = front{2+i};
        end
        if strcmp(T4files{i},'')
            T4files{i} = front{3+numsubjects+i};
        end
    end
    
    for subject = 1:length(glmfiles)
        cd(outputdir)
        subnameindex = strfind(glmfiles{subject},'/vc');
        
        slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
        
        subname{subject} = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
        
        [glmparamname glmparamval] = textread(glmfiles{subject},'%s%s',20,'delimiter',':');
        
        concfile = glmparamval{strcmp('name of data file',deblank(glmparamname))};
        slashloc = strfind(concfile,'/');
        concfolder = concfile(slashloc(1):slashloc(end));
        if ~strcmp(concfile(slashloc(1):slashloc(1)+4),'/data')
            glmslashloc = strfind(glmfiles{subject},'/');
            concfolder = glmfiles{subject}(1:glmslashloc(end));
        end
        
        if ~isempty(strfind(list,'28Subjects_RePreprocessed'))
            concfolder = ['/data/cn3/joe/ResourceDataLimited/' subname{subject} '/'];
        end
        
        concfile = concfile(slashloc(end)+1:end);
        
        if stevedata
            concfolder(1:length('/data/nil-external/mcl/Nelson')) = [];
            concfolder = ['/data/cn3/steven' concfolder];
        end
        
        badlocation = strfind(concfolder,'Create_GLMv1R/..');
        if ~isempty(badlocation)
            concfolder = [concfolder(1:badlocation-1) concfolder(badlocation + length('Create_GLMv1R/..') : end)];
        end
        
        [concfront concruns] = textread([concfolder concfile],'%s%s',20,'delimiter',':');
        concruns = concruns(strcmp('file',deblank(concfront)));
        FDformattot = [];
        
        newconcfilename = [outputdir '/goodvoxels/' subname{subject} '.conc'];
        delete([newconcfilename]);
        fid = fopen([newconcfilename],'at'); %open the output file for writing
        fclose(fid);
        dlmwrite([newconcfilename],['number_of_files:   ' num2str(length(concruns))],'-append','delimiter','');
        
        for run = 1:length(concruns)
            if stevedata
                concruns{run}(1:length('/data/nil-external/mcl/Nelson')) = [];
                concruns{run} = ['/data/cn3/steven' concruns{run}];
            end    
            
                slashloc = strfind(concruns{run},'/');
                subnameloc = strfind(concruns{run},subname{subject});
                movement = [concruns{run}(1:(subnameloc(1)+length(subname{subject}))) 'movement/' concruns{run}(slashloc(end)+1:end-9) '.dat'];
                %movement = [concfolder 'movement/' concruns{run}(slashloc(end)+1:end-9) '.dat'];
            
            
            
            [drmstotal drmstrans drmsrot drmscol ddt_mvm] = rdat_calculations(movement,4,50);
            delta_mvm = [0 0 0 0 0 0 ; (ddt_mvm(2:end,:) - ddt_mvm(1:end-1,:))];
            FD=(sum(abs(delta_mvm),2));
            FDformat = repmat('x',1,length(FD));
            FDformat(logical(FD<.2)) = '+';
            FDformat(1:4) = 'x';
            FDformattot = [FDformattot FDformat];
            
            dlmwrite([newconcfilename],['file:' concruns{run}(1:end-9) '_333.4dfp.img'],'-append','delimiter','');
            
        end
        
        dlmwrite('./goodvoxels/total_tmask_avi.txt',FDformattot,'')
        
        clear FDformattot
        
        
        
        %if ~exist(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subname{subject} '/7112b_fs_LR/Ribbon/'])
        if ~exist([outputdir '/' subname{subject} '/7112b_fs_LR/Ribbon/'])
            mkdir([outputdir '/' subname{subject} '/7112b_fs_LR/Ribbon/'])
            system(['csh /data/cn4/evan/Scripts/create_ribbon_task.csh ' subname{subject} ' ' outputdir])
        end
        
        system(['csh /data/cn4/evan/Scripts/RibbonVolumetoSurfaceMapping_taskdata.csh ' subname{subject}])
        copyfile([outputdir 'goodvoxels/' subname{subject} '_goodvoxels.nii.gz'],[outputdir 'goodvoxels/temp' subname{subject} '_goodvoxels.nii.gz'])
        delete([outputdir 'goodvoxels/' subname{subject} '*'])
        copyfile([outputdir 'goodvoxels/temp' subname{subject} '_goodvoxels.nii.gz'],[outputdir 'goodvoxels/' subname{subject} '_goodvoxels.nii.gz'])
        delete([outputdir 'goodvoxels/temp*'])
        
        surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR';
        HEMS = {'L';'R'};
        smooth = '2.55';
        
        conditionfiles = dir([directory subname{subject} '*.4dfp.img']);
        if jessdata
            conditionfiles = [conditionfiles ; dir(['/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/' subname{subject} '*.4dfp.img'])];
        end
            
        if isempty(conditionfiles)
            
            %subname{subject} = ['_' glmfiles{subject}];
            undsubname{subject} = ['_' glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1)];
            conditionfiles = dir([directory undsubname{subject} '*.4dfp.img']);
        end
        
        submask = [outputdir 'goodvoxels/' subname{subject} '_goodvoxels.nii.gz'];
        
%         for condition = 1:length(conditionfiles)
%             
%             subfunc = [directory '/' conditionfiles(condition).name(1:end-9)];
%             if jessdata && (condition > length(dir([directory subname{subject} '*.4dfp.img'])))
%                 subfunc = ['/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/' conditionfiles(condition).name(1:end-9)];
%             end
%             system(['niftigz_4dfp -n ' subfunc ' ' subfunc])
%             
%             for hem = 1:2
%                 
%                 midsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.midthickness.native.surf.gii'];
%                 midsurf_LR32k = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
%                 whitesurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.white.native.surf.gii'];
%                 pialsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.pial.native.surf.gii'];
%                 nativedefsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
%                 outsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
%         
%                 surfname = [subname{subject} '_' HEMS{hem} '_condition' num2str(condition)];
%                 disp('Map volume to surface')
%                 system([workbenchdir '/wb_command -volume-to-surface-mapping ' subfunc '.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
%                 disp('Dilate surface timecourse')
%                 system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
%                 
%                 disp('Deform timecourse to 32k fs_LR')
%                 cd([surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/'])
%                 %system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' outputdir '/' surfname '_dil10.func.gii ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii']);
%                 system([workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
%                 
%                 disp('Smooth surface timecourse')
%                 system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii ' smooth ' ' outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii'])
%                 
%                 system(['rm ' outputdir '/' surfname '.func.gii'])
%                 system(['rm ' outputdir '/' surfname '_dil10.func.gii'])
%                 system(['rm ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii'])
%             end
%         end
        
    end
    clear subname
end
