
cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list';
outputdir = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/voxelwiseFC/';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
ROImetric = 'minima_4mm_rois.L.func.gii';
ROIdir = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/';
BOLD_suffix = '_BOLD_L_smooth2.55_32k_fsLR';

Hem = 'L';
output_singlesubs = 0;

%Cohortfile format:
% subjectname funcpath/funcvol
%
% e.g.
%
% vc33416 /data/cn4/laumannt/vc33416_rest1.4dfp.img


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

roistruct = gifti_evan([ROIdir '/' ROImetric]);


% Make output folder
system(['mkdir ' outputdir])



for s = 1:subnum
    tmask = load(tmaskfiles{s});
    
        BOLDname = [subjects{s} BOLD_suffix];
        system(['caret_command64 -file-convert -format-convert ASCII ' BOLDname '.func.gii']);
        system(['rm ' BOLDname '_noHEAD.func.gii']);
        system(['awk ''NF > 25'' ' BOLDname '.func.gii > ' BOLDname '_noHEAD.func.gii']);
        
        % Load surface and volume data
        surf_BOLD = load([BOLDname '_noHEAD.func.gii']);
        surf_BOLD(:,1) = [];
        system(['rm ' BOLDname '_noHEAD.func.gii']);
        surf_BOLD = surf_BOLD(:,logical(tmask));

    
    for roinum = 1:length(roistruct.data)
        
        roiname{roinum} = roistruct.data{roinum}.metadata(1).value;
        
        disp(['Subject ' subjects{s} ', ' roiname{roinum}])
        
        seedtimecourse = mean(surf_BOLD(logical(roistruct.data{roinum}.data),:),1);
        
        fcmap{roinum}(:,s) = FisherTransform(corr(seedtimecourse',surf_BOLD'));
        
        if output_singlesubs
            
            save(gifti(fcmap{roinum}(:,s)),[outputdir '/' subjects{s} '_' roiname{roinum} '_' Hem '.func.gii'],'ExternalFileBinary');
                        
        end
        
    end
    
end

for roinum = 1:length(roistruct.data)
    
    avgfcmap = mean(fcmap{roinum},2);
    save(gifti(avgfcmap),[outputdir '/Avg_' roiname{roinum} '_' Hem '.func.gii'],'ExternalFileBinary');
    
end









