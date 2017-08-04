outputdir = '/data/cn4/evan/Task_parcellation/Manytasks_subregister/';


lists = {'/data/cn4/evan/Task_parcellation/SteveSubjects/28subs_16tps_Ret_list_of_glms.glm_list',...
    '/data/cn4/evan/Task_parcellation/SteveSubjects/28subs_16tps_Prim_list_of_glms.glm_list',...
    '/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list',...
    '/data/cn4/evan/Occipitotemporal/Orthography/Adults.glm_list',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
    '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt',...
    '/data/cn3/joe/SustainedTaskLoad/FloatSustainedOnly_Scrubbed_f06_28Subjects.list',...
    '/data/cn4/maital/Errors/2TaskGroups/ACRN.glm_list',...
    '/data/cn4/maital/reward/REPREPROCESS/33Ss_NotScrubbed/2x3x5_list_of_glms.glm_list'};

directories = {'/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/SlowReveal_Oldnew/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/SlowReveal_Priming/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/Katie/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/Orthography/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/VisualSwitch/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/ResourceDataLimited/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/VisualAttention/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/SustainedTaskLoad/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/ACRN/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/Reward/SingleSubjects/'};%,...

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';

hems = {'L' 'R'};



for listnum = 4%1:length(lists)
    list = lists{listnum};
    disp(list)
    directory = directories{listnum};
    
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
        disp(glmfiles{subject})
        cd(outputdir)
        subnameindex = strfind(glmfiles{subject},'/vc');
        
        slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
        
        subname = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
        
        surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR';
        HEMS = {'L';'R'};
        smooth = '2.55';
        
        conditionfiles = dir([directory subname '*.4dfp.img']);
        if isempty(conditionfiles)
            conditionfiles = dir([directory '_' subname '*.4dfp.img']);
        end
        
        submask = [outputdir 'goodvoxels/' subname '_goodvoxels.nii.gz'];
        %submask = ['/data/cn4/evan/Task_parcellation/Manytasks_subregister/goodvoxels/' subname{subject} '_goodvoxels.nii.gz'];
        
        for condition = 1:length(conditionfiles)
            
            subfunc = [directory '/' conditionfiles(condition).name(1:end-9)];
            system(['niftigz_4dfp -n ' subfunc ' ' subfunc])
            
            if ~exist([subname '.nii.gz'])
                copyfile([subfunc '.nii.gz'],[subname '.nii.gz'])
            else
                system(['fslmerge -t ' subname '.nii.gz ' subname '.nii.gz ' subfunc '.nii.gz'])
            end
        end
        
        gunzip([subname '.nii.gz']);
        
        for hem = 1:2
            
            surfname = [subname '_' HEMS{hem}];
            
            %if ~exist([outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii']) %|| ~isempty(strfind(list,'28subs_16tps_Ret_list_of_glms')) || ~isempty(strfind(list,'28subs_16tps_Prim_list_of_glms'))
            
            midsurf = [surfdir '/' subname '/7112b_fs_LR/Native/' subname '.' HEMS{hem} '.midthickness.native.surf.gii'];
            midsurf_LR32k = [surfdir '/' subname '/7112b_fs_LR/fsaverage_LR32k/' subname '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
            whitesurf = [surfdir '/' subname '/7112b_fs_LR/Native/' subname '.' HEMS{hem} '.white.native.surf.gii'];
            pialsurf = [surfdir '/' subname '/7112b_fs_LR/Native/' subname '.' HEMS{hem} '.pial.native.surf.gii'];
            nativedefsphere = [surfdir '/' subname '/7112b_fs_LR/Native/' subname '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
            outsphere = [surfdir '/' subname '/7112b_fs_LR/fsaverage_LR32k/' subname '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
            
            disp('Map volume to surface')
            system([workbenchdir '/wb_command -volume-to-surface-mapping ' subname '.nii ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
            disp('Dilate surface timecourse')
            system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
            
            disp('Deform timecourse to 32k fs_LR')
            cd([surfdir '/' subname '/7112b_fs_LR/fsaverage_LR32k/'])
            %system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' outputdir '/' surfname '_dil10.func.gii ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii']);
            system([workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
            
            disp('Smooth surface timecourse')
            system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii ' smooth ' ' outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii'])
            
            %end
            
            %thisdata = gifti([outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii']);
            %Mergedata{hem} = [Mergedata{hem} thisdata.cdata];
            
            system(['rm ' outputdir '/' surfname '.func.gii'])
            system(['rm ' outputdir '/' surfname '_dil10.func.gii'])
            system(['rm ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii'])
            
            cd(outputdir)
        end
        
        
        
        system(['wb_command -cifti-create-dense-timeseries ' subname '.dtseries.nii -volume ' subname '.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' subname '_L_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/hcp-zfs/home/laumannt/120_parcellation/L.atlasroi_group_noproj.func.gii -right-metric ' subname '_R_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/hcp-zfs/home/laumannt/120_parcellation/R.atlasroi_group_noproj.func.gii'])
        
        delete([outputdir '/' subname '.nii'])
        
        delete([outputdir '/' subname '_L_dil10_32k_fsLR_smooth' smooth '.func.gii'])
        delete([outputdir '/' subname '_R_dil10_32k_fsLR_smooth' smooth '.func.gii'])
        clear subname
    end
    
    %system(fslmergestr)
    
end
delete([outputdir '/vc*.nii.gz'])





