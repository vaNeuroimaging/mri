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
hemname_low = {'left';'right'};

numsubs = 0;

Mergedata{1} = [];
Mergedata{2} = [];
for listnum = 1:length(lists)
    
    list = lists{listnum};
    
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
        
        numsubs = numsubs+1;
        
        fslmergestr = 'fslmerge -t Temp.nii.gz';
        disp(['List number ' num2str(listnum) ', Subject number ' num2str(subject) ', overall subject number ' num2str(numsubs)])
        cd(outputdir)
        subnameindex = strfind(glmfiles{subject},'/vc');
        
        slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
        
        subname{subject} = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
        
        surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR';
        HEMS = {'L';'R'};
        smooth = '2.55';
        
        %         conditionfiles = dir([directory subname{subject} '*.4dfp.img']);
        %         if isempty(conditionfiles)
        %             conditionfiles = dir([directory '_' subname{subject} '*.4dfp.img']);
        %         end
        %
        %         submask = [outputdir 'goodvoxels/' subname{subject} '_goodvoxels.nii.gz'];
        %         %submask = ['/data/cn4/evan/Task_parcellation/Manytasks_subregister/goodvoxels/' subname{subject} '_goodvoxels.nii.gz'];
        %
        %         for condition = 1:length(conditionfiles)
        %             fslmergestr = [fslmergestr ' ' directory conditionfiles(condition).name(1:end-9) '.nii.gz'];
        %         end
        %         system(fslmergestr)
        %
        %         gunzip([outputdir '/Temp.nii.gz'])
        %
        %         for hem = 1:2
        %
        %             %surfname = [subname{subject} '_' HEMS{hem}];
        %             surfname = ['Temp_' HEMS{hem}];
        %
        %
        %             midsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.midthickness.native.surf.gii'];
        %             midsurf_LR32k = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        %             whitesurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.white.native.surf.gii'];
        %             pialsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.pial.native.surf.gii'];
        %             nativedefsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        %             outsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        %
        %             disp('Map volume to surface')
        %             system([workbenchdir '/wb_command -volume-to-surface-mapping ' outputdir '/Temp.nii ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
        %             disp('Dilate surface timecourse')
        %             system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
        %
        %             disp('Deform timecourse to 32k fs_LR')
        %             cd([surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/'])
        %             %system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' outputdir '/' surfname '_dil10.func.gii ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii']);
        %             system([workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
        %
        %             disp('Smooth surface timecourse')
        %             system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii ' smooth ' ' outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii'])
        %
        %         end
        %
        %         disp('Creating dense timeseries')
        %         if ~exist([outputdir '/' subname{subject} '.dtseries.nii'])
        %             system(['wb_command -cifti-create-dense-timeseries ' outputdir '/' subname{subject} '.dtseries.nii -volume ' outputdir '/Temp.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' outputdir '/Temp_L_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii -right-metric ' outputdir '/Temp_R_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii'])
        %         else
        %             system(['wb_command -cifti-create-dense-timeseries ' outputdir '/' subname{subject} '2.dtseries.nii -volume ' outputdir '/Temp.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric ' outputdir '/Temp_L_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii -right-metric ' outputdir '/Temp_R_dil10_32k_fsLR_smooth' smooth '.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii'])
        %             system(['wb_command -cifti-merge ' outputdir '/Temp.dtseries.nii -cifti ' outputdir '/' subname{subject} '.dtseries.nii -cifti ' outputdir '/' subname{subject} '2.dtseries.nii'])
        %             system(['mv ' outputdir '/Temp.dtseries.nii ' outputdir '/' subname{subject} '.dtseries.nii'])
        %             delete([outputdir '/' subname{subject} '2.dtseries.nii'])
        %         end
        
        
        
        for hem = 1:2
            disp(['Conducting correlation for ' HEMS{hem} ' hemisphere'])
            system([workbenchdir '/wb_command -cifti-correlation ' outputdir '/' subname{subject} '.dtseries.nii ' outputdir '/Temp_corr.dconn.nii -roi-override -' hemname_low{hem} '-roi /data/cn4/laumannt/subcortical_mask/' HEMS{hem} '.atlasroi_erode3.32k_fs_LR.shape.gii -fisher-z'])
            
            system([workbenchdir '/wb_command -cifti-math "a" ' outputdir '/Temp2_corr.dconn.nii -fixnan 0 -var a ' outputdir '/Temp_corr.dconn.nii'])
            
            system(['mv ' outputdir '/Temp2_corr.dconn.nii ' outputdir '/Temp_corr.dconn.nii'])
            
            if subject==1%numsubs==1%(listnum==1) && (s == 1)
                disp('Initializing average correlation matrix')
                system(['mv ' outputdir '/Temp_corr.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
            else
                disp('Adding to running average')
                system([workbenchdir '/wb_command -cifti-math "a+b" ' outputdir '/avg_temp.dconn.nii -var a ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii -var b ' outputdir '/Temp_corr.dconn.nii']);
                disp('done')
                system(['mv ' outputdir '/avg_temp.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
                system(['rm ' outputdir '/Temp_corr.dconn.nii'])
            end
        end
        delete([outputdir '/Temp*']);
    end
    
    
    for hem = 1:2
        disp('Calculating average correlation map')
        system([workbenchdir '/wb_command -cifti-math "a/' num2str(length(glmfiles)) '" ' outputdir '/avg_temp.dconn.nii -var a ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii '])
        disp('done')
        system(['mv ' outputdir '/avg_temp.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
        
        
        if listnum==1
            disp('Initializing task average correlation matrix')
            system(['mv ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii' outputdir '/taskavg_corr_' HEMS{hem} '.dconn.nii'])
        else
            disp('Adding to running task average')
            system([workbenchdir '/wb_command -cifti-math "a+b" ' outputdir '/avg_temp.dconn.nii -var a ' outputdir '/taskavg_corr_' HEMS{hem} '.dconn.nii -var b ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii']);
            disp('done')
            system(['mv ' outputdir '/avg_temp.dconn.nii ' outputdir '/taskavg_corr_' HEMS{hem} '.dconn.nii'])
            system(['rm ' outputdir '/avg_corr.dconn.nii'])
        end
        
    end
    
    
    
    
end



