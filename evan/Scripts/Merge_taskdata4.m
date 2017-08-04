clear
warning off

outname = '/data/cn4/evan/Task_parcellation/Manytasks_subregister/Manytasks';

lists = {'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/VisualAttention_fbfdc05_f0_b0.list',...
    '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};%,...

%'/data/cn4/evan/Task_parcellation/SteveSubjects/SteveFakeList_newdata.glm_list',...};
    

directories = {'/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Trials/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/Script/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/Joe3Tasks_IndConditions/SingleSubjects/'};%,...

%'/data/cn4/evan/Task_parcellation/SteveSubjects/OrigData/',...};
    
%
smoothnum = 2.55;
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
Hems = {'L' 'R'};

extranamestuffstrings = {'_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_'};

%[],...};
    

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
        
        try
            
            subnameindex = strfind(glmfiles{subject},'/vc');
            
            slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
            
            subname{subject} = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
            
            conditionfiles = dir([directory subname{subject} '*.4dfp.img']);
            
        catch
            
            subname{subject} = [glmfiles{subject}];
            conditionfiles = dir([directory '_' subname{subject} '*.4dfp.img']);
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
        
        condition_surfdataL = [];
        condition_surfdataR = [];
        condition_voldata = [];
        
        subjectcounter = 0;
        
        for subject = 1:length(glmfiles)
            
            file = dir([directory subname{subject} allconditionnames{conditionnum}]);
            if isempty(file)
                file = dir([directory '_' subname{subject} allconditionnames{conditionnum}]);
            end
            
            if ~isempty(file)
                
                disp(subname{subject})
                
                surfdir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subname{subject} '/7112b_fs_LR/'];
                
                if exist(surfdir)
                    
                    subjectcounter = subjectcounter+1;
                    
                    condition_voldata(:,1:timepoints(conditionnum),subjectcounter) = read_4dfpimg([directory file(1).name]);
                    
                    system(['nifti_4dfp -n ' directory file(1).name ' ' outname 'Temp.nii']);
                    
                    for hemnum = 1:2;
                        
                        midsurf = [surfdir '/Native/' subname{subject} '.' Hems{hemnum} '.midthickness.native.surf.gii'];
                        whitesurf = [surfdir '/Native/' subname{subject} '.' Hems{hemnum} '.white.native.surf.gii'];
                        pialsurf = [surfdir '/Native/' subname{subject} '.' Hems{hemnum} '.pial.native.surf.gii'];
                        
                        system([workbenchdir '/wb_command -volume-to-surface-mapping ' outname 'Temp.nii ' midsurf ' ' outname 'Temp' Hems{hemnum} '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -voxel-subdiv 5']);
                        system([workbenchdir '/wb_command -metric-smoothing ' midsurf ' ' outname 'Temp' Hems{hemnum} '.func.gii ' num2str(smoothnum) ' ' outname 'Temp_smooth_' Hems{hemnum} '.func.gii'])
                        
                        
                        % Deform from native to 32K fs_LR surface
                        cd([surfdir '/fsaverage_LR32k/'])
                        system(['caret_command64 -deformation-map-apply native232k_fs_LR.' Hems{hemnum} '.deform_map METRIC_AVERAGE_TILE '  outname 'Temp_smooth_' Hems{hemnum} '.func.gii '  outname 'Temp_smooth_' Hems{hemnum} '_32k_fsLR.func.gii']);
                        
                        system(['caret_command64 -file-convert -format-convert XML ' outname 'Temp_smooth_' Hems{hemnum} '_32k_fsLR.func.gii'])
                        %system(['awk ''NF > 25'' ' outname 'Temp_smooth_' Hems{hemnum} '_32k_fsLR.func.gii > ' outname 'Temp_smooth_' Hems{hemnum} '_32k_fsLR_noHEAD.func.gii'])
                        
                    end
                    
                    tempdata = gifti([outname 'Temp_smooth_L_32k_fsLR.func.gii']);
                    condition_surfdataL(:,1:timepoints(conditionnum),subject) = tempdata.cdata;
                    
                    tempdata = gifti([outname 'Temp_smooth_R_32k_fsLR.func.gii']);
                    condition_surfdataR(:,1:timepoints(conditionnum),subject) = tempdata.cdata;
                    
                    delete([outname 'Temp*']);
                    
                else
                    disp(['Subject ' subname{subject} ' has no fsLR. Skipping'])
                end
            end
            
            
        end
        
        voldata = [voldata mean(condition_voldata,3)];
        surfdataL = [surfdataL mean(condition_surfdataL,3)];
        surfdataR = [surfdataR mean(condition_surfdataR,3)];
        
        clear condition_voldata condition_surfdataL condition_surfdataR
        
    end
    
    clear conditionname
end

write_4dfpifh([outname '.4dfp.ifh'], size(voldata,2),'bigendian');
write_4dfpimg(voldata,[outname '.4dfp.img'],'bigendian');
save(gifti(single(surfdataL)),[outname '_L.func.gii']);
save(gifti(single(surfdataR)),[outname '_R.func.gii']);

