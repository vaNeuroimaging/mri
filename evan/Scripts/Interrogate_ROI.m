clear
warning off

%list = '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list';
%'/data/cn4/evan/Occipitotemporal/VisualAttention/VisualAttention_fbfdc05_f0_b0.list';
%'/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list';
%
%
%timepoints = 9;
%
%directory = '/data/cn4/evan/Occipitotemporal/VisualSwitch/Script/SingleSubjects/';
%'/data/cn4/evan/Occipitotemporal/VisualAttention/Script/SingleSubjects/';
%'/data/cn4/evan/Test/SingleSubjects/';
%
%

lists = {'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CplusT_SCRUBfbf07_all30_10tps_020312.glm_list'};
%{'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
 %   '/data/cn4/evan/Occipitotemporal/VisualAttention/VisualAttention_fbfdc05_f0_b0.list',...
  %  '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list'};

directories = {'/data/cn4/evan/Occipitotemporal/VisualSwitch/Script_CplusT/SingleSubjects/'};
%{'/data/cn4/evan/Occipitotemporal/VisualSwitch/Script/SingleSubjects/',...
%    '/data/cn4/evan/Occipitotemporal/VisualAttention/Script/SingleSubjects/',...
%    '/data/cn4/evan/Test/SingleSubjects/'};

%listtimepoints = [9 7 7];


ROIs = {'lOccip_TypexTime_ROI'};
    
roidir = '/data/cn4/evan/Occipitotemporal/VisualSwitch/Script_CplusT/';
%'/data/cn4/evan/Occipitotemporal/VisualSwitch/Script/';
% roinames = dir([roidir '/ROI*.nii']);
% for i = 1:length(roinames)
%     ROIs{i} = roinames(i).name(1:end-4);
% end



outputfilename = [roidir '/Timecourses.txt'];

extranamestuffstrings = {'_avgtc_'...
    '_avgtc_'...
    '_avgtc_'};


delete([outputfilename]);
fid = fopen([outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Condition','ROI','Timepoint','Value'); %write the output file header
fclose(fid);
dlmwrite([outputfilename],' ','-append');

for roinum = 1:length(ROIs);
    try roidata{roinum} = read_4dfpimg_HCP([roidir ROIs{roinum} '.4dfp.img']);
    catch
        evalc(['!nifti_4dfp -4 ' roidir ROIs{roinum} '.nii ' roidir ROIs{roinum} '.4dfp.img']);
        roidata{roinum} = read_4dfpimg_HCP([roidir ROIs{roinum} '.4dfp.img']);
    end
%     thisdata = load_nii([directory ROIs{roinum}]);
%     roidata{roinum} = reshape(thisdata.img,[size(thisdata.img,1) * size(thisdata.img,2) * size(thisdata.img,3) , 1]);
end


for listnum = 1:length(lists)
list = lists{listnum};
directory = directories{listnum};
%timepoints = listtimepoints(listnum);
extranamestuff = extranamestuffstrings{listnum};


[front entries] = textread(list,'%s%s','delimiter',':');

numsubjects = str2num(entries{2});
glmfiles = entries(3:2+numsubjects);
T4files = entries(4+numsubjects:end);
for i = 1:length(glmfiles)
    if strcmp(glmfiles{i},'')
        glmfiles{i} = front{2+i};
    end
    if strcmp(T4files{i},'')
        T4files{i} = front{3+numsubjects+i};
    end
end





for subject = 1:length(glmfiles)
    %for timepoint = 1:timepoints
    
    subnameindex = strfind(glmfiles{subject},'/vc');
    
    slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
    
    subname = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
    
        
        conditionfiles = dir([directory subname '*.4dfp.img']);
        %conditionfiles = dir([directory  'avgtc_*.4dfp.img']);
        
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
                
                conditionname = conditionfiles(condition).name((length(subname)+length(extranamestuff) + 1):conditionnameend);
            
                disp(['Subject ' subname ', ' conditionname])
                
            
            imagedata = read_4dfpimg_HCP([directory conditionfiles(condition).name]);
            
            
            
            for roinum = 1:length(ROIs)
                
%                 %----------------------------------------------------------
%                 roimask = roidata{roinum};
%                 roiindices = find(roimask);
%                 roimask(:) = 0;
%                 roimask(roiindices(randi(length(roiindices)))) = 1;
%                 %----------------------------------------------------------
                
                for timepoint = 1:timepoints
                
                maskedmean = mean(imagedata(roidata{roinum}>0,timepoint));
                
                %maskedmean = mean(imagedata(roimask>0,timepoint));
                
                texttowrite = [subname,'   ',conditionname,'   ',ROIs{roinum},'   ',num2str(timepoint),'   ',num2str(maskedmean)];  %save the data as a string to be written to the output
                
                dlmwrite([outputfilename],texttowrite,'-append','delimiter','');%write the data to the output file
                end
                
            end
            
        end
    %end
end
end