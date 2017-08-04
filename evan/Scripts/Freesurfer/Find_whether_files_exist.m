clear
warning off


lists = {'/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list'};
%     '/data/cn4/evan/Occipitotemporal/Orthography/Adults.glm_list',...
%     '/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt',...
% '/data/cn4/evan/Task_parcellation/SteveSubjects/Subjectlist_withanat.txt',...
% '/data/cn4/maital/reward/REPREPROCESS/33Ss_NotScrubbed/2x3x5_list_of_glms.glm_list',...
% '/data/cn4/maital/Errors/2TaskGroups/ACRN.glm_list',...
% '/data/cn3/joe/SustainedTaskLoad/FloatSustainedOnly_Scrubbed_f06_28Subjects.list',...
% '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
% '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
% '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};%,...
    


directories = {'/data/cn4/evan/Occipitotemporal/Orthography/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualAttention/SingleSubjects/',...
    '/data/cn4/evan/Task_parcellation/SteveSubjects/OrigData/',...
    '/data/cn4/evan/Occipitotemporal/Reward/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/ACRN/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/SustainedTaskLoad/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Trials/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/',...
    '/data/cn4/evan/Occipitotemporal/Joe3Tasks_IndConditions/SingleSubjects/'};%,...
    %


Hems = {'L' 'R'};

extranamestuffstrings = {'_avgtc_' ...
    '_avgtc_' ...
    [],...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_' ...
    '_avgtc_'};

outputfilename = '/data/cn4/evan/Scripts/Freesurfer/freesurfer_datalist.txt';
delete([outputfilename]);
fid = fopen([outputfilename],'at'); %open the output file for writing
fclose(fid);
dlmwrite([outputfilename],' ','-append');


for listnum = 1:length(lists)
    
    if listnum == 4
        1;
    end
    
    list = lists{listnum};
    disp(list)
    %directory = directories{listnum};
    %timepoints = listtimepoints(listnum);
    %extranamestuff = extranamestuffstrings{listnum};
    
    [front entries] = textread(list,'%s%s','delimiter',':');
    
    numsubjects = str2num(entries{2});
    glmfiles = entries(3:2+numsubjects);
    T4files = entries(4+numsubjects:end);
    for i = 1:length(glmfiles)
        if strcmp(glmfiles{i},'')
            glmfiles{i} = front{2+i};
        end
        %     if strcmp(T4files{i},'')
        %         T4files{i} = front{3+numsubjects+i};
        %     end
    end
    
    for i = 1:length(T4files)
        if strcmp(T4files{i},'')
            T4files{i} = front{3+numsubjects+i};
        end
    end
    
        
        for subject = 1:length(glmfiles)
            

            
                subnameindex = strfind(glmfiles{subject},'/vc');
            
                slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
            
                subname = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
                

            
            
            %if exist(['/data/cn4/segmentation/freesurfer5_supercomputer/' subname '/' subname '_mpr_n1_111_t88.nii.gz']);
            if  exist(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subname '/7112b_fs_LR/fsaverage_LR32k/' subname '.R.midthickness.32k_fs_LR.surf.gii']);
                %
                %
            %if exist(['/data/cn4/segmentation/' subname '/']) | exist(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subname '/7112b_fs_LR/']);
                disp([subname ': yes'])
            else
                disp([subname ': no'])
                
                 subnameindex = strfind(T4files{subject},'/vc');
                 anatfolder = T4files{subject}(1:subnameindex-1);
%                 disp(T4files{subject})
%                 
                 dlmwrite([outputfilename],[anatfolder '   ' subname],'-append','delimiter','');%write the data to the output file
                
            end
        end
end
            
            
%             
%             
%             
%             conditionfiles = dir([directory subname '*.4dfp.img']);
%             
%             
%             for condition = 1:length(conditionfiles)
%                 
%                 timestr = [];
%                 for t=1:20
%                     timestr = [timestr '_' num2str(t) ];
%                     
%                     if isempty(strfind(conditionfiles(condition).name,timestr))
%                         break
%                     end
%                     
%                     conditionnameend = strfind(conditionfiles(condition).name,timestr) - 1;
%                     timepoints = t;
%                     
%                 end
%                 
%                 %conditionnameend = strfind(conditionfiles(condition).name,timestr) - 1;
%                 
%                 conditionname = conditionfiles(condition).name((length(subname)+length(extranamestuff) + 1):conditionnameend);
%                 
%                 if ~exist([directory conditionfiles(condition).name(1:end-9) '_' Hem '.func.gii'])
%                     
%                     map_vol_to_surface([directory conditionfiles(condition).name],Hem);
%                     
%                 end
%                 
%                 imagedata = gifti([directory conditionfiles(condition).name(1:end-9) '_' Hem '.func.gii']);
%                 imagedata = imagedata.cdata;
%                 
%                 
%                 
%                 for roinum = 1:size(roidata,2)
%                     
%                     if size(roidata,2)>1
%                         roiname{roinum} = [ROImetrics{roifilenum}(1:end-9) '_map' num2str(roinum)];
%                     else
%                         roiname{roinum} = [ROImetrics{roifilenum}(1:end-9)];
%                     end
%                     
%                     
%                     disp(['Subject ' subname ', ' conditionname ', ' roiname{roinum}])
%                     
%                     
%                     for timepoint = 1:timepoints
%                         
%                         datainROI = imagedata(logical(roidata(:,roinum)),timepoint);
%                         
%                         datainROI = datainROI(~logical(isnan(datainROI)));
%                         
%                         maskedmean = mean(datainROI);
%                         
%                         texttowrite = [subname,'   ',conditionname,'   ',roiname{roinum},'   ',num2str(timepoint),'   ',num2str(maskedmean)];  %save the data as a string to be written to the output
%                         
%                         dlmwrite([outputfilename],texttowrite,'-append','delimiter','');%write the data to the output file
%                         
%                     end
%                     
%                     
%                     
%                 end
%                 
%             end
%             
%         end
%     end
% end
