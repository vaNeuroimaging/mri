
cohortfile = '/data/cn4/evan/RestingState/FC_Mapping/28Subjects_surfaceparcellation.list';
outputdir = '/data/cn4/evan/RestingState/FC_Mapping/lOT/lOT_minima_fc';
tmasklist = '/data/cn4/evan/RestingState/28subjects_FDp2_00_DVinit_20_00/28subjects_FDp2_00_DVinit_20_00_newDV_tmasklist.txt';
seeds = {'DA_ant','DA_post','DMN_temp','FPC','Vis'};
seedfolder = '/data/cn4/evan/RestingState/FC_Mapping/lOT/lOT_minima_fc/';
output_singlesubs = 0;

%Cohortfile format:
% subjectname funcpath/funcvol
%
% e.g.
%
% vc33416 /data/cn4/laumannt/vc33416_rest1.4dfp.img

mask4dfp='/data/cn4/evan/ROIs/glm_atlas_mask_333.4dfp.img';
etype = 'bigendian';

[ign subnum] = system(['cat ' cohortfile ' | wc -l']);
subnum = str2num(subnum);

% Read in subject names, functional volume locations
for s = 1:subnum
    [ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
    [ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
    funcpaths{s} = funcpaths{s}(1:end-1);
    [ign funcvols{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1}''']);
    funcvols{s} = funcvols{s}(1:end-1);
end

[tmasksubjects tmaskfiles]=textread(tmasklist,'%s%s');
if ~isequal(tmasksubjects,subjects)
    error('masklist subjects do not match cohortfile subjects');
end
for i=1:numel(tmaskfiles)
    fprintf('\t%d\t%s\n',i,subjects{i,1});
    if ~exist(tmaskfiles{i,1})
        error('\t%s is not found\n',tmaskfiles{i,1});
    end
end

% Make output folder
system(['mkdir ' outputdir])
mask = read_4dfpimg(mask4dfp);

for seed = 1:length(seeds)
    try seeddata{seed} = read_4dfpimg([seedfolder seeds{seed} '.4dfp.img']);
    catch
        evalc(['!nifti_4dfp -4 ' seedfolder seeds{seed} '.nii ' seedfolder seeds{seed} '.4dfp.img']);
        seeddata{seed} = read_4dfpimg([seedfolder seeds{seed} '.4dfp.img']);
    end
    
end

for s = 1:subnum
    tmask = load(tmaskfiles{s});
    BOLD = read_4dfpimg([funcpaths{s} '/' funcvols{s} '.4dfp.img']);
    BOLD = BOLD(:,logical(tmask));
    
    for seed = 1:length(seeds)
        
        disp(['Subject ' subjects{s} ', ' seeds{seed}])
        
        seedtimecourse = mean(BOLD(logical(seeddata{seed}),:),1);
        
        fcmap{seed}(:,s) = FisherTransform(corr(seedtimecourse',BOLD'));
        fcmap{seed}(~logical(mask),s) = 0;
        
        if output_singlesubs
            
            write_4dfpimg(fcmap{seed}(:,s),fullfile(outputdir,[subjects{s} '_' seeds{seed} '.4dfp.img']),etype);
            write_4dfpifh(fullfile(outputdir,[subjects{s} '_' seeds{seed} '.4dfp.img']),1,etype);
            
        end
        
    end
    
end

for seed = 1:length(seeds)
    avgfcmap = mean(fcmap{seed},2);
    write_4dfpimg(avgfcmap,fullfile(outputdir,['Avg_' seeds{seed} '.4dfp.img']),etype);
    write_4dfpifh(fullfile(outputdir,['Avg_' seeds{seed} '.4dfp.img']),1,etype);
end









