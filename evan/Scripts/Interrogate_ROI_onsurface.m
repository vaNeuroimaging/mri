clear
warning off

%list = '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list';
%'/data/cn4/evan/Occipitotemporal/VisualAttention/VisualAttention_fbfdc05_f0_b0.list';
%'/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list';
%
%
%

lists = {'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list'};
%{'/data/cn4/evan/Occipitotemporal/Lexicalization/n28_D1_glm_event_cln_gen_anal_wt4.list'};
%{'/data/cn4/evan/Task_parcellation/SteveSubjects/SteveFakeList_newdata.glm_list'};
%{'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
 %   '/data/cn4/evan/Occipitotemporal/VisualAttention/VisualAttention_fbfdc05_f0_b0.list',...
  %  '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list',...
   % '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};%,...
    %


directories = {'/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Cues/SingleSubjects/'};
%{'/data/cn4/evan/Occipitotemporal/Lexicalization/SingleSubjects/'};
%{'/data/cn4/evan/Task_parcellation/SteveSubjects/OrigData/'};
% {'/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Trials/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/VisualAttention/Script/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/GlassPattern/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/MRT/SingleSubjects/'};%,...
    %


ROImetrics = {'ConsensusLH_FPCDA_noedge_clusters_1_to_100_minsize40.func.gii','ConsensusRH_FPCDA_noedge_clusters_1_to_100_minsize40.func.gii'};
ROIdir = '/data/cn4/evan/RestingState/Consensus/';

Hems = {'L' 'R'};




outputfilename = ['/data/cn4/evan/Occipitotemporal/Consensus_Edgeconstrained_Timecourses.txt'];

extranamestuffstrings = {[]};
% {'_avgtc_'...
%     '_avgtc_' ...
%     '_avgtc_' ...
%     '_avgtc_'};



delete([outputfilename]);
fid = fopen([outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Condition','ROI','Timepoint','Value'); %write the output file header
fclose(fid);
dlmwrite([outputfilename],' ','-append');


for listnum = 1:length(lists)
    list = lists{listnum};
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
    
    for roifilenum = 1:length(ROImetrics)
        
        Hem = Hems{roifilenum};
        
        roidata = gifti([ROIdir '/' ROImetrics{roifilenum}]);
        roidata = roidata.cdata;
        
        
        
        for subject = 1:length(glmfiles)
            
            try
            
                subnameindex = strfind(glmfiles{subject},'/vc');
            
                slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
            
                subname = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
                
            catch
            
                subname = ['_' glmfiles{subject}];
            end
            
            
            conditionfiles = dir([directory subname '*.4dfp.img']);
            
            
            for condition = 1:length(conditionfiles)
                
                timestr = [];
                for t=1:20
                    timestr = [timestr '_' num2str(t) ];
                    
                    if isempty(strfind(conditionfiles(condition).name,timestr))
                        break
                    end
                    
                    conditionnameend = strfind(conditionfiles(condition).name,timestr) - 1;
                    timepoints = t;
                    
                end
                
                %conditionnameend = strfind(conditionfiles(condition).name,timestr) - 1;
                
                conditionname = conditionfiles(condition).name((length(subname)+length(extranamestuff) + 1):conditionnameend);
                
                if ~exist([directory conditionfiles(condition).name(1:end-9) '_' Hem '.func.gii'])
                    
                    map_vol_to_surface([directory conditionfiles(condition).name],Hem);
                    
                end
                
                imagedata = gifti([directory conditionfiles(condition).name(1:end-9) '_' Hem '.func.gii']);
                imagedata = imagedata.cdata;
                
                
                
                for roinum = 1:size(roidata,2)
                    
                    if size(roidata,2)>1
                        roiname{roinum} = [ROImetrics{roifilenum}(1:end-9) '_map' num2str(roinum)];
                    else
                        roiname{roinum} = [ROImetrics{roifilenum}(1:end-9)];
                    end
                    
                    
                    disp(['Subject ' subname ', ' conditionname ', ' roiname{roinum}])
                    
                    
                    for timepoint = 1:timepoints
                        
                        datainROI = imagedata(logical(roidata(:,roinum)),timepoint);
                        
                        datainROI = datainROI(~logical(isnan(datainROI)));
                        
                        maskedmean = mean(datainROI);
                        
                        texttowrite = [subname,'   ',conditionname,'   ',roiname{roinum},'   ',num2str(timepoint),'   ',num2str(maskedmean)];  %save the data as a string to be written to the output
                        
                        dlmwrite([outputfilename],texttowrite,'-append','delimiter','');%write the data to the output file
                        
                    end
                    
                    
                    
                end
                
            end
            
        end
    end
end
