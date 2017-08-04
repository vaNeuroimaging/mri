clear
warning off

outname = '/data/cn4/evan/Task_parcellation/Manytasks/Manytasks';

lists = {%'/data/cn4/evan/Occipitotemporal/Lexicalization/n28_D1_glm_event_cln_gen_anal_wt4.list',...
    %'/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list',...
    '/data/cn4/evan/Occipitotemporal/Orthography/Adults.glm_list',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt',...
'/data/cn3/steven/NP847/glmsLists/28subs_16tps_Prim_list_of_glms.glm_list',...'/data/cn4/evan/Task_parcellation/SteveSubjects/Subjectlist_withanat.txt',...
'/data/cn4/maital/reward/REPREPROCESS/33Ss_NotScrubbed/2x3x5_list_of_glms.glm_list',...
'/data/cn4/maital/Errors/2TaskGroups/ACRN.glm_list',...
'/data/cn3/joe/SustainedTaskLoad/FloatSustainedOnly_Scrubbed_f06_28Subjects.list',...
'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
'/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};

directories = {%'/data/cn4/evan/Occipitotemporal/Lexicalization/SingleSubjects/',...
    %'/data/cn4/evan/Occipitotemporal/Katie/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/Orthography/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/SteveSubjects/OrigData/',...
    '/data/cn4/evan/Occipitotemporal/Reward/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/ACRN/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/SustainedTaskLoad/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Trials/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/Joe3Tasks_IndConditions/SingleSubjects/'};%,...

extranamestuffstrings = {%'_avgtc_' ...
    %'_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    [],...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_'};


Hems = {'L' 'R'};


voldata = [];
surfdataL = [];
surfdataR = [];


for listnum = 1:length(lists)
    list = lists{listnum};
    disp(list)
    directory = directories{listnum};
    %timepoints = listtimepoints(listnum);
    extranamestuff = extranamestuffstrings{listnum};
    
    [front entries] = textread(list,'%s%s','delimiter',':');
    
    numsubjects = str2num(entries{2});
    glmfiles = entries(3:2+numsubjects);
    %T4files = entries(4+numsubjects:end);
    for i = 1:length(glmfiles)
        if strcmp(glmfiles{i},'')
            glmfiles{i} = front{2+i};
        end
        %     if strcmp(T4files{i},'')
        %         T4files{i} = front{3+numsubjects+i};
        %     end
    end
    
    allconditionnames = {[' ']};
    
    for subject = 1:length(glmfiles)
        
        
            
            subnameindex = strfind(glmfiles{subject},'/vc');
            
            slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
            
            subname{subject} = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
            
            conditionfiles = dir([directory subname{subject} '*.4dfp.img']);
            
        if isempty(conditionfiles)
            
            %subname{subject} = ['_' glmfiles{subject}];
            subname{subject} = ['_' glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1)];
            conditionfiles = dir([directory subname{subject} '*.4dfp.img']);
        end
        
        
        
        
        
        for condition = 1:length(conditionfiles)
            
            timestr = [];
            for t=1:20
                timestr = [timestr '_' num2str(t) ];
                
                if isempty(strfind(conditionfiles(condition).name,timestr))
                    break
                end
                
                conditionnameend = strfind(conditionfiles(condition).name,timestr) - 1;
                timepoints(condition) = t;
                
            end
            
            %conditionnameend = strfind(conditionfiles(condition).name,timestr) - 1;
            
            %conditionname{condition,subject} = [directory conditionfiles(condition).name((length(subname)+length(extranamestuff) + 1):conditionnameend)];
            
            thiscondition = [conditionfiles(condition).name((length(subname{subject}) + 1):end)];
            
            match = 0;
            for i = 1:length(allconditionnames)
                if strcmp(allconditionnames{i},thiscondition)
                    match = 1;
                end
            end
            if match==0
                allconditionnames{end+1} = thiscondition;
            end
            
        end
    end
    allconditionnames = allconditionnames(2:end);
    
    
    for conditionnum = 1:length(allconditionnames)
        disp(['Condition ' allconditionnames{conditionnum}])
        
        %condition_surfdataL = [];
        %condition_surfdataR = [];
        condition_voldata = [];
        
        for subject = 1:length(glmfiles)
            
            file = dir([directory subname{subject} allconditionnames{conditionnum}]);
            
            if ~isempty(file)
            
            
            disp(subname{subject})
            
            condition_voldata(:,1:timepoints(conditionnum),subject) = read_4dfpimg([directory file(1).name]);
            
%             if ~exist([directory file(1).name(1:end-9) '_L.func.gii'])
%                 map_vol_to_surface([directory file(1).name],'L');
%             end
%             tempdata = gifti([directory file(1).name(1:end-9) '_L.func.gii']);
%             condition_surfdataL(:,1:timepoints(conditionnum),subject) = tempdata.cdata;
%             
%             if ~exist([directory file(1).name(1:end-9) '_R.func.gii'])
%                 map_vol_to_surface([directory file(1).name],'R');
%                 delete([directory file(1).name(1:end-9) '.nii'])
%             end
%             tempdata = gifti([directory file(1).name(1:end-9) '_R.func.gii']);
%             condition_surfdataR(:,1:timepoints(conditionnum),subject) = tempdata.cdata;
            
            end
            
        end
        
        voldata = [voldata mean(condition_voldata,3)];
        %surfdataL = [surfdataL mean(condition_surfdataL,3)];
        %surfdataR = [surfdataR mean(condition_surfdataR,3)];
        
        clear condition_voldata condition_surfdataL condition_surfdataR
        
    end
    
    clear conditionname
end

write_4dfpifh([outname '.4dfp.ifh'], size(voldata,2),'bigendian');
write_4dfpimg(voldata,[outname '.4dfp.img'],'bigendian');
% save(gifti(single(surfdataL)),[outname '_L.func.gii']);
% save(gifti(single(surfdataR)),[outname '_R.func.gii']);

