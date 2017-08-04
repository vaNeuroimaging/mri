outputdir = '/data/cn4/evan/Task_parcellation/Manytasks_subregister/';
cd(outputdir)

lists = {%'/data/cn4/evan/Occipitotemporal/Lexicalization/n28_D1_glm_event_cln_gen_anal_wt4.list',...
    %'/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list',...
    %'/data/cn4/evan/Occipitotemporal/Orthography/Adults.glm_list',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt',...
    '/data/cn4/evan/Task_parcellation/SteveSubjects/Subjectlist_withanat.txt',...'/data/cn3/steven/NP847/glmsLists/28subs_16tps_Prim_list_of_glms.glm_list',...
    '/data/cn4/maital/reward/REPREPROCESS/33Ss_NotScrubbed/2x3x5_list_of_glms.glm_list',...
    '/data/cn4/maital/Errors/2TaskGroups/ACRN.glm_list',...
    '/data/cn3/joe/SustainedTaskLoad/FloatSustainedOnly_Scrubbed_f06_28Subjects.list',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...    '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
    '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};

subject = 0;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = ~mask.cdata;

for listnum = 5:length(lists)
    list = lists{listnum};
    disp(list)
    
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
        
        
        subnameindex = strfind(glmfiles{subject},'/vc');
        
        slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
        
        subname{subject} = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
        
        subfiles = dir([subname{subject} '*.func.gii']);
        data = gifti(subfiles(1).name); data = data.cdata;
        
        numzeros = nnz(abs(data(logical(mask),:))<.00001);
        numnans = nnz(isnan(data(logical(mask),:)));
        
        if numzeros
            disp([num2str(numzeros) ' zeros in subject ' subname{subject}])
        end
        if numnans
            disp([num2str(numnans) ' NaNs in subject ' subname{subject}])
        end
        
    end
    clear subname
end


