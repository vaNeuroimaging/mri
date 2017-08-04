outputdir = '/data/cn4/evan/Task_parcellation/Manytasks_subregister/';


lists = {%'/data/cn4/evan/Occipitotemporal/Lexicalization/n28_D1_glm_event_cln_gen_anal_wt4.list',...
    %'/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list',...
    %'/data/cn4/evan/Occipitotemporal/Orthography/Adults.glm_list',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt',...
    '/data/cn4/evan/Task_parcellation/SteveSubjects/Subjectlist_withanat2.txt',...'/data/cn3/steven/NP847/glmsLists/28subs_16tps_Prim_list_of_glms.glm_list',...
    '/data/cn4/maital/reward/REPREPROCESS/33Ss_NotScrubbed/2x3x5_list_of_glms.glm_list',...
    '/data/cn4/maital/Errors/2TaskGroups/ACRN.glm_list',...
    '/data/cn3/joe/SustainedTaskLoad/FloatSustainedOnly_Scrubbed_f06_28Subjects.list',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...    '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
    '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};

directories = {%'/data/cn4/evan/Occipitotemporal/Lexicalization/SingleSubjects/',...
    %'/data/cn4/evan/Occipitotemporal/Katie/SingleSubjects/',...
    %'/data/cn4/evan/Occipitotemporal/Orthography/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/SteveSubjects/OrigData/',...
    '/data/cn4/evan/Occipitotemporal/Reward/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/ACRN/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/SustainedTaskLoad/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Trials/SingleSubjects/',...    '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/Joe3Tasks_IndConditions/SingleSubjects/'};%,...

extranamestuffstrings = {%'_avgtc_' ...
    %'_avgtc_' ...
    %'_avgtc_' ...
    '_avgtc_' ...
    [],...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...    '_avgtc_' ...
    '_avgtc_'};

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';

hems = {'L' 'R'};


Mergedata{1} = [];
Mergedata{2} = [];
for listnum = 1:length(lists)
    fslmergestr = 'fslmerge -t Manytasks.nii.gz Manytasks.nii.gz';
    list = lists{listnum};
    disp(list)
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
        
        conditionfiles = dir([directory subname{subject} '*.nii.gz']);
        if jessdata
            conditionfiles = [conditionfiles ; dir(['/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/' subname{subject} '*.4dfp.img'])];
        end
        if isempty(conditionfiles)
            
            %subname{subject} = ['_' glmfiles{subject}];
            undsubname{subject} = ['_' glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1)];
            conditionfiles = dir([directory undsubname{subject} '*.nii.gz']);
        end
        
        for hem = 1:2
            for condition = 1:length(conditionfiles)
                
                if hem==1
                    
                    if listnum == 1 && subject == 1 && condition == 1
                        copyfile([directory conditionfiles(condition).name],'Manytasks.nii.gz')
                    else
                        
                        if jessdata && (condition > length(dir([directory subname{subject} '*.4dfp.img'])))
                            fslmergestr = [fslmergestr ' /data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/' conditionfiles(condition).name(1:end-9)];
                            %system(['fslmerge -t Manytasks.nii.gz Manytasks.nii.gz /data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/' conditionfiles(condition).name(1:end-9)]);
                        else
                            %system(['fslmerge -t Manytasks.nii.gz Manytasks.nii.gz ' directory conditionfiles(condition).name])
                            fslmergestr = [fslmergestr ' ' directory conditionfiles(condition).name];
                        end
                    end
                end
            
 %               conditiondata = gifti([outputdir subname{subject} '_' hems{hem} '_condition' num2str(condition) '_dil10_32k_fsLR_smooth2.55.func.gii']);
 %               Mergedata{hem} = [Mergedata{hem} conditiondata.cdata];
                
            end
        end
        
    end
    clear subname
    system(fslmergestr)
end

%save(gifti(single(Mergedata{1})),'Manytasks_L.func.gii')
%save(gifti(single(Mergedata{2})),'Manytasks_R.func.gii')
%
gunzip Manytasks.nii.gz
system('wb_command -cifti-create-dense-timeseries Manytasks.dtseries.nii -volume Manytasks.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric Manytasks_L.func.gii -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii -right-metric Manytasks_R.func.gii -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii')
