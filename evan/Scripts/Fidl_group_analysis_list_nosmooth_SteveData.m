

%User enter values here
%=========================================================================

%GLM List
list = '/data/cn4/evan/Task_parcellation/SteveSubjects/28subs_16tps_Ret_list_of_glms.glm_list';
%'/data/cn4/evan/Task_parcellation/SteveSubjects/28subs_16tps_Prim_list_of_glms.glm_list';
%'/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list';
%'/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list';
%'/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt';
%'/data/cn3/joe/SustainedTaskLoad/FloatSustainedOnly_Scrubbed_f06_28Subjects.list';
%'/data/cn4/maital/Errors/2TaskGroups/ACRN.glm_list';
%'/data/cn4/maital/reward/REPREPROCESS/33Ss_NotScrubbed/2x3x5_list_of_glms.glm_list';
%
%'/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list';
%'/data/cn4/evan/Occipitotemporal/Orthography/Adults.glm_list';
%'/data/cn4/evan/Occipitotemporal/Lexicalization/n28_D1_glm_event_cln_gen_anal_wt4.list';
%;
%
%
%
%
%
%
%
%

%Name of factors (other than time) in anova design. If time is the only
%factor, enter 'Time'.
factornames = {'Temp'};

%ANOVAsetup reflects the design of the anova. Each factor is represented as
%the dimension of this variable that factor varies along. So factor 1 above
%represents the differences between the rows of this variable, while factor
%2 represents the differences between the columns.

namestofind = {'HIT4' 'CR4' 'HIT5' 'CR5' 'HIT6' 'CR6' 'HIT7' 'CR7'};
conditionstorun = {[1 2] [3 4] [5 6] [7 8]};

%ANOVAsetup = {1 2 3 4 5 6 7 8 9 10 11 12};

%Katie?
%ANOVAsetup = {[1 3 5 7] [9]; [2 4 6 8] [10]};

%Reward
% ANOVAsetup(:,:,1) = {[1 2] [6 7]; [11 12] [16 17]; [21 22] [26 27]};
% %ANOVAsetup(:,:,2) = {[2] [7]; [12] [17]; [22] [27]};
% %ANOVAsetup(:,:,3) = {[3] [8]; [13] [18]; [23] [28]};
% ANOVAsetup(:,:,2) = {[4 5] [9 10]; [14 15] [19 20]; [24 25] [29 30]};
% %ANOVAsetup(:,:,5) = {[5] [10]; [15] [20]; [25] [30]};

%;


%levelnames is a variable of the same outer cell dimensionality as factornames
%and the same inner cell dimensionality as ANOVAsetup. It contains cells
%with multiple character strings representing the names of the levels
%within each factor. 
levelnames = {'Temp'};


%The directory that results will be written to.
outputdir = '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/SlowReveal_Oldnew/';

%Do you want to also generate single-subject timecourses for each cell of
%the anova? (hint: this saves later work)
RunSingleSubs = 1;

%Space that the output will be written into
voxelspace = '333';

%Prefix written out before output files.
outputname = 'ANOVA';

%Thresholds for monte-carlo correction
montecarloZ = 3.75;
montecarlok = 18;

%=========================================================================



warning off


%Determine linux version
temp = computer;
if strcmp(temp(end-1:end),'64')
    linuxdir = '/home/usr/fidl/fidl_code/fidl_2.64/bin_linux64';
else
    linuxdir = '/home/usr/fidl/fidl_code/fidl_2.64/bin_linux';
end

%Specify directory (inside the outputdir) that temp files will write to
Scratchdir = 'SCRATCH/';

%Read the GLM list file
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

%Make the output directories
mkdir(outputdir);
mkdir([outputdir '/' Scratchdir]);
if RunSingleSubs
    mkdir([outputdir '/SingleSubjects/']);
end

%Create headers for the driver file
driverfile{1} = ['subject'];
for factor = 1:length(factornames)
    driverfile{1} = [driverfile{1}  '    ' factornames{factor}];
end
if numel(conditionstorun) > 1
    driverfile{1} = [driverfile{1} '    time'];
end
driverfile{1} = [driverfile{1}  '    *.4dfp.img'];


%Create the driver file and write out the headers
driverfilename = [outputdir outputname '_driver.dat'];
delete(driverfilename);
fid = fopen(driverfilename,'at'); %open the output file for writing
fprintf(fid,'%s\n',driverfile{1}); %write the output file header
fclose(fid);
 
%Start making the command that will run the ANOVA
anovarunstr = ['!' linuxdir '/fidl_anova_new -driver "' driverfilename '" -voxel_threshold 0.01 -uncompress /home/usr/fidl/lib/glm_atlas_mask_' voxelspace '.4dfp.img -output Z_uncorrected Z_monte_carlo  -glmfiles '];


for subject = 1:length(glmfiles)
    
    clear conditionindex
    
    %Figure out the subject name from the glm file name (as long as it starts with vc)
    subnameindex = strfind(glmfiles{subject},'/vc');
    slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
    subname = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
    
    %Figure out if the glm file name has anything else in addition to the subject name
    glmsuffix = glmfiles{subject}(subnameindex(end)+1+length(subname):end-4);
    
    
    disp(['Subject ' subname])
    
    
    
    
    
    for length = 500:-1:10
        fid = fopen(glmfiles{subject},'r+');
     try
        header = textscan(fid,'%s',length,'Delimiter','%\n');
        fclose(fid);
        break
     catch
         fclose(fid);
     end
    end
    clear length
    
     header = header{1};
     conditionnum = 0;
     for line = 1:size(header,1)
         if ~isempty(strfind(header{line},'glm effect label'))
             conditionnum = conditionnum+1;
         for namenum = 1:size(namestofind,2)
         
            if ~isempty(strfind(header{line},namestofind{namenum}))
                
                namefound(namenum) = conditionnum;
            end
         end
         end
     end
     
     clear ANOVAsetup
     
     for icondition = 1:size(conditionstorun,2)
         ANOVAsetup{icondition} = [];
         for isublevel = 1:length(conditionstorun{icondition})
             if (length(namefound) >= conditionstorun{icondition}(isublevel)) && (namefound(conditionstorun{icondition}(isublevel))>0)
                ANOVAsetup{icondition} = [ANOVAsetup{icondition} namefound(conditionstorun{icondition}(isublevel))];
             end
         end
         
     end
     
     condstokeep = ones(1,length(ANOVAsetup));
     for icondition = 1:length(ANOVAsetup)
         if isempty(ANOVAsetup{icondition})
             condstokeep(icondition) = 0;
         end
     end
     ANOVAsetup = ANOVAsetup(logical(condstokeep));
         
    
     clear length header conditionnum namefound
    
    
    
    %This section reads the glm file header and figures out the names,
    %design matrix column numbers, and durations of the conditions
    %---------------------------------------------------------------------
    conditionnamelines = evalc(['!grep ''glm effect label'' ' glmfiles{subject} ' -a']);
    conditionnameindex1 = strfind(conditionnamelines,':= ');
    conditionnameindex2 = strfind(conditionnamelines,'glm effect label');
    
    conditionlengthlines = evalc(['!grep ''glm effect length'' ' glmfiles{subject} ' -a']);
    conditionlengthindex1 = strfind(conditionlengthlines,':= ');
    conditionlengthindex2 = strfind(conditionlengthlines,'glm effect length');
    
    conditionindexlines = evalc(['!grep ''glm effect column'' ' glmfiles{subject} ' -a']);
    conditionindexindex1 = strfind(conditionindexlines,':= ');
    conditionindexindex2 = strfind(conditionindexlines,'glm effect column');
    
    for condition = 1:length(conditionnameindex1)-1
        conditionnames{condition} = conditionnamelines(conditionnameindex1(condition)+3:conditionnameindex2(condition+1)-2);
        conditionlengths(condition) = str2num(conditionlengthlines(conditionlengthindex1(condition)+3:conditionlengthindex2(condition+1)-2));
        conditionindex(condition) = str2num(conditionindexlines(conditionindexindex1(condition)+3:conditionindexindex2(condition+1)-2));
        
    end
    conditionnames{length(conditionnameindex1)} = conditionnamelines(conditionnameindex1(end)+3:end-1);
    conditionlengths(length(conditionlengthindex1)) = str2num(conditionlengthlines(conditionlengthindex1(end)+3:end-1));
    conditionindex(length(conditionindexindex1)) = str2num(conditionindexlines(conditionindexindex1(end)+3:end-1));
    
%     if subject == 1
%         sub1conditionnames = conditionnames;
%         sub1conditionlengths = conditionlengths;
%         sub1conditionindex = conditionindex;
%     else
%         
%         baseconditionnames = conditionnames;
%         baseconditionlengths = conditionlengths;
%         baseconditionindex = conditionindex;
%         
%         for condition = 1:length(baseconditionnames)
%             newindex = find(strcmp(baseconditionnames{condition},sub1conditionnames));
%             if length(newindex)==1 && ~(newindex == condition)
%                 conditionnames{newindex} = baseconditionnames{condition};
%                 conditionlengths(newindex) = baseconditionlengths(condition);
%                 conditionindex(newindex) = baseconditionindex(condition);
%             end
%         end
%     end
                
    
    conditionindex = conditionindex +1;
    %---------------------------------------------------------------------
        
    
    %Start making the command that will run the single-subject aspect of
    %the overall anova
    subjzstatrunstr = ['!' linuxdir '/fidl_zstat -glm_file ' glmfiles{subject} ' -tc '];
    
    %Start making the command that will create single-subject timecourses
    singlesubjrunstr = ['!' linuxdir '/fidl_avg_zstat -glm_files ' glmfiles{subject} ' -tc '];
    
    for cell = 1:numel(ANOVAsetup)
        timelengths(cell) = min(conditionlengths(ANOVAsetup{cell}));
    end
    timepoints = min(timelengths);
    timepoints = repmat(timepoints,[1,numel(ANOVAsetup)]);
    
    
    %This section loops through conditions/timepoints and adds pertinent information
    %about them to the driver file and to the relevant commands
    %---------------------------------------------------------------------
     for cell = 1:numel(ANOVAsetup)
         string = '[';
         
         %Figure out the number of factors
         if isempty(find(size(ANOVAsetup)==1))
            nfactors = ndims(ANOVAsetup);
         else
             nfactors = 1;
         end
         
         %For each factor, add a column to the driver file
         for factor = 1:nfactors;
             string = [string 'dim{' num2str(factor) '} '];
         end
         string = [string '] = ind2sub(size(ANOVAsetup),cell)'];
         evalc(string);
         
         %Start figuring out the name of the temp file that will be written to the
         %scratch directory for this subject/factor
         singlesubjectoutfilename{cell} = ['"' subname '_avgtc_' conditionnames{ANOVAsetup{cell}}];

         %Figure out the number of timepoint needed for this factor's
         %timecourse
         %timepoints(cell) = max(conditionlengths(ANOVAsetup{cell}));
         
         
            for time = 1:timepoints%(cell)
                
                %add timepoint numbers to the temp file name
                singlesubjectoutfilename{cell} = [singlesubjectoutfilename{cell} '_' num2str(time)];
                
                %Figure out the column numbers for each timepoint
                timecourse{subject}{time}{cell} = [];
                
                %Add subject number to the driver file
                driverfile{end+1} = num2str(subject); 
                
                %Add level number to the driver file
                if numel(ANOVAsetup)>1
                    string = [];
                    for factor = 1:length(dim)
                        
                %        driverfile{end} = [driverfile{end} '    ' levelnames{factor}{dim{factor}}];
                        
                        %string = [string num2str(dim{factor}) ','];
                    end
                    %string = string(1:end-1);
                    %evalc(['driverfile{end} = [driverfile{end} ''    '' levelnames{' string '}]']);
                end
                
                %Add timepoint to the driver file
                driverfile{end} = [driverfile{end}  '    ' num2str(time)];
                
                %Add the temp file name to the driver file
                driverfile{end} = [driverfile{end} '    ' Scratchdir subname glmsuffix '_'];
                for condnum = 1:length(ANOVAsetup{cell})
                    driverfile{end} = [driverfile{end} conditionnames{ANOVAsetup{cell}(condnum)} '+'];
                    timecourse{subject}{time}{cell} = [timecourse{subject}{time}{cell} num2str(conditionindex(ANOVAsetup{cell}(condnum)) +time -1) '+'];
                end
                driverfile{end} = [driverfile{end}(1:end-1) '_' num2str(time) '.4dfp.img'];
                                
                %Make the timecourse for this subject/condition
                timecourse{subject}{time}{cell} = [timecourse{subject}{time}{cell}(1:end-1) ' '];
                
                %Add relevent timecourse numbers to the commands that will
                %be run
                subjzstatrunstr = [subjzstatrunstr timecourse{subject}{time}{cell} ' '];
                singlesubjrunstr = [singlesubjrunstr timecourse{subject}{time}{cell} ' '];
                
                %Write out the driver file
                dlmwrite(driverfilename,driverfile{end},'-append','delimiter','');
                
            end
            
            
            singlesubjectoutfilename{cell} = [singlesubjectoutfilename{cell} '"'];
                
            
     end
     %---------------------------------------------------------------------
     
     %Finish making the single subject ANOVA and timecourse commands
     subjzstatrunstr = [subjzstatrunstr '-xform_file ' T4files{subject} ' -compress /home/usr/fidl/lib/glm_atlas_mask_' voxelspace '.4dfp.img -atlas ' voxelspace ' -scratchdir ' outputdir '/' Scratchdir];
     singlesubjrunstr = [singlesubjrunstr '-xform_files ' T4files{subject} ' -mask /home/usr/fidl/lib/glm_atlas_mask_' voxelspace '.4dfp.img -atlas ' voxelspace ' -directory ' outputdir '/SingleSubjects/  -frames ' num2str(timepoints) ' -glmpersub 1 -tc_names '];
     
     %Add the temp filenames to the single subject ANOVA command
     for cell = 1:numel(ANOVAsetup)
         singlesubjrunstr = [singlesubjrunstr singlesubjectoutfilename{cell} ' '];
     end
     
     %Run the single subject ANOVA command
     sublog{subject} = evalc(subjzstatrunstr);
     
     %Run the single subject timecourse command
     if RunSingleSubs
         singlesublog{subject} = evalc(singlesubjrunstr);
     end
     
     %Add this subject's glm filename to the overall ANOVA command
     anovarunstr = [anovarunstr glmfiles{subject} ' '];
     
     clear conditionnames conditionlengths conditionindex baseconditionnames baseconditionlengths baseconditionindex conditionnameindex1 conditionnameindex2 conditionlengthindex1 conditionlengthindex2 conditionindexindex1 conditionindexindex2
    
end

% %Add the relevent timecourses for each subject/condition/timepoint to the overall ANOVA command
% anovarunstr = [anovarunstr '-tc '];
% for cell = 1:numel(ANOVAsetup)
% for time = 1:timepoints%(cell)
%     for subject = 1:length(glmfiles)
%         anovarunstr = [anovarunstr timecourse{subject}{time}{cell}(1:end-1) ','];
%     end
%     anovarunstr = [anovarunstr(1:end-1) ' '];
% end
% end
% 
% %Change to the output directory
% cd(outputdir)
% 
% disp('Group ANOVA')
% 
% %Finish making the overall ANOVA command
% anovarunstr = [anovarunstr '-Nimage_name "' outputdir outputname '_Nimage.4dfp.img" -threshold_extent "' num2str(montecarloZ) ' ' num2str(montecarlok) '" -pval .05 -scratchdir ' outputdir '/' Scratchdir ' -GIGAdesign -glmpersub ' num2str(ones(1,length(glmfiles))) ' -clean_up'];
% 
% %Run the overall ANOVA command
% anovalog = evalc(anovarunstr);
%      
% 
% outputfiles = dir([outputdir '/*.4dfp.img']);
% for filenum=1:length(outputfiles)
%     if ~exist([outputdir '/' outputfiles(filenum).name(1:end-9) '.nii']);
%         evalc(['!nifti_4dfp -n ' outputdir '/' outputfiles(filenum).name ' ' outputdir '/' outputfiles(filenum).name(1:end-9) '.nii']);
%     end
% end


    