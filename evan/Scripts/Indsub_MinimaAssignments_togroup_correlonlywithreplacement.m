minimadist = 3;
costparameter = 2; %relative weighting of distance vs connectivity pattern

hem = 'L';
cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list';
edgemapfolder = '/data/cn4/laumannt/left_hem/';
groupedgemapfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/MinimaAssignments/';
groupedgefile = 'avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg.func.gii';
tmaskname = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
medialmask = '/data/cn4/evan/RestingState/FC_Mapping_120/medialmask.func.gii';

makenewgroupstuff = 0;

outputstatsfile = [outputfolder '/Statistics.txt'];

delete(outputstatsfile);
fid = fopen(outputstatsfile,'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','GroupNode','SubjectNode','CostType','Cost'); %write the output file header
fclose(fid);
dlmwrite(outputstatsfile,' ','-append');

surfacecoordfile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.coord.gii';
surfacetopofile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.32k_fs_LR.topo.gii';

medialmaskdata = gifti(medialmask);
medialmaskind = find(~medialmaskdata.cdata);

groupminimafile = [groupedgefile(1:end-9) '_minima_dist' num2str(minimadist) '.func.gii'];

disp('Calculating group minima')
if makenewgroupstuff
    metric_minima([groupedgemapfolder '/' groupedgefile],minimadist,[],hem,outputfolder,groupminimafile);
end

groupminimadata = gifti([outputfolder '/' groupminimafile]);
groupminimadata = groupminimadata.cdata;
groupminimadata(medialmaskind) = 0;
groupminimavertices = find(groupminimadata);

disp('Calculating distances from group minima')
if makenewgroupstuff
    
    distances = zeros(length(groupminimavertices),length(groupminimadata));
    
    for minimanum = 1:length(groupminimavertices)
        
        string{minimanum} = ['    Minimum number ' num2str(minimanum) ' of ' num2str(length(groupminimavertices))];
        if minimanum==1; fprintf('%s',string{minimanum}); else fprintf([repmat('\b',1,length(string{minimanum-1})) '%s'],string{minimanum}); end
        
        system(['caret_command64 -surface-geodesic ' surfacecoordfile ' ' surfacetopofile ' temp.func.gii false -node ' num2str(groupminimavertices(minimanum)-1)]);
        evalc('!caret_command64 -file-convert -format-convert XML temp.func.gii');
        thisdistances = gifti('temp.func.gii');
        distances(minimanum,:) = thisdistances.cdata;
    end
    disp(' ')
    save([outputfolder '/distances'],'distances')
    system('rm temp.func.gii');
else
    load([outputfolder '/distances']);
end

maxdist = max(max(distances));

[tmasksubjects tmaskfiles]=textread(tmaskname,'%s%s');

[ign subnum] = system(['cat ' cohortfile ' | wc -l']);
subnum = str2num(subnum);

avgcorrelmat = zeros(length(groupminimavertices),length(groupminimadata));


disp('Calculating average connectivity patterns of group minima');

for s = 1:subnum
    [ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
    [ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
    [ign funcvols{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1}''']);
    [ign surfdirs{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' cohortfile]);
    
    surfBOLDfilename = dir([groupedgemapfolder '/' subjects{s} '_BOLD*.func.gii']);
    surfBOLDfile{s} = [groupedgemapfolder '/' surfBOLDfilename(1).name(1:end-9)];
    
    tmask{s} = load(tmaskfiles{s});
    
    if makenewgroupstuff
        disp(['    Subject ' subjects{s}])
        evalc(['!rm ' surfBOLDfile{s} '_noHEAD.func.gii']);
        evalc(['!awk ''NF > 25'' ' surfBOLDfile{s} '.func.gii > ' surfBOLDfile{s} '_noHEAD.func.gii']);
        surf_BOLD = load([surfBOLDfile{s} '_noHEAD.func.gii']);
        surf_BOLD(:,1) = [];
        surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
        evalc(['!rm ' surfBOLDfile{s} '_noHEAD.func.gii']);
        
        correlmat = FisherTransform(paircorr_mod(surf_BOLD(groupminimavertices,:)',surf_BOLD'));
        correlmat(:,medialmaskind) = 0;
        
        avgcorrelmat = avgcorrelmat + correlmat;
    end
    
end

if makenewgroupstuff
    avgcorrelmat = avgcorrelmat ./ subnum;
    save([outputfolder '/avgcorrelmat'],'avgcorrelmat');
else
    load([outputfolder '/avgcorrelmat']);
end




disp('Conducting minima assignments');
for s = 1:subnum
    
    disp(['Subject ' subjects{s} ':']);
    disp('    calculating minima');
    subjectedgemetric = dir([edgemapfolder '/' subjects{s} '_avg_edge_avg_smooth_' hem '_noalone.func.gii']);
    subjectedgemetric = subjectedgemetric(1).name;
    
    subjectminimametric = [subjectedgemetric(1:end-9) '_minima_dist' num2str(minimadist) '.func.gii'];
    
    try subjectminimadata = gifti([outputfolder '/' subjectminimametric]);
    catch
        metric_minima([edgemapfolder '/' subjectedgemetric],minimadist,[],hem,outputfolder,subjectminimametric);
        subjectminimadata = gifti([outputfolder '/' subjectminimametric]);
    end
    subjectminimadata = subjectminimadata.cdata;
    subjectminimadata(medialmaskind) = 0;
    subjectminimavertices = find(subjectminimadata);
    
    disp('    calculating connectivity of minima');
    
    evalc(['!awk ''NF > 25'' ' surfBOLDfile{s} '.func.gii > ' surfBOLDfile{s} '_noHEAD.func.gii']);
    surf_BOLD = load([surfBOLDfile{s} '_noHEAD.func.gii']);
    surf_BOLD(:,1) = [];
    surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
    evalc(['!rm ' surfBOLDfile{s} '_noHEAD.func.gii']);
    
    subjcorrelmat = FisherTransform(paircorr_mod(surf_BOLD(subjectminimavertices,:)',surf_BOLD'));
    subjcorrelmat(:,medialmaskind) = 0;
    
    
    disp('    calculating assignments');
    
    subvsgroupcorrel = paircorr_mod(avgcorrelmat',subjcorrelmat');
    
    [ign assignment] = max(subvsgroupcorrel,[],2);
    
     subjvsgroupdist = distances(:,subjectminimavertices);
%     
%     subjvsgroupcost = (subjvsgroupdist ./ maxdist) .* costparameter  + (1-(subvsgroupcorrel+1)/2);
%     
%     [assignment, cost] = munkres(subjvsgroupcost);
    
    for minimanum = 1:length(groupminimavertices)
        
        if s==1; outputmetric = zeros(size(groupminimadata),length(groupminimavertices)); end
        outputmetric(subjectminimavertices(assignment(minimanum)),minimanum) = 1;
        outputmetric(groupminimavertices(minimanum),minimanum) = 30;
        
        texttowrite = [subjects{s} '   ' num2str(groupminimavertices(minimanum)) '   ' num2str(subjectminimavertices(assignment(minimanum))) '   Distance   '   num2str(subjvsgroupdist(minimanum,assignment(minimanum)))];
        dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
        
        texttowrite = [subjects{s} '   ' num2str(groupminimavertices(minimanum)) '   ' num2str(subjectminimavertices(assignment(minimanum))) '   Correl   '   num2str(subvsgroupcorrel(minimanum,assignment(minimanum)))];
        dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
        
%         texttowrite = [subjects{s} '   ' num2str(groupminimavertices(minimanum)) '   ' num2str(subjectminimavertices(assignment(minimanum))) '   TotalCost_costparam' num2str(costparameter) '   ' num2str(subjvsgroupcost(minimanum,assignment(minimanum)))];
%         dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
    end
    
end

disp('Writing output');

save(gifti(outputmetric),[outputfolder '/SubjectsMinimaAssignments_dist' num2str(minimadist) '_correlonlywithreplacement.func.gii']);

% outputstring = ['wb_command -metric-merge ' outputfolder '/SubjectsMinimaAssignments_dist' num2str(minimadist) '_param' num2str(costparameter) '.func.gii'];
% 
% for minimanum = 1:length(groupminimavertices)
%     save(gifti(outputmetric{minimanum}),[outputfolder '/TempSubjectsMinimaAssignments_' num2str(minimanum) '.func.gii'],'ExternalFileBinary');
%     outputstring = [outputstring ' -metric ' outputfolder '/TempSubjectsMinimaAssignments_' num2str(minimanum) '.func.gii'];
% end
% system(outputstring);
% system(['rm ' outputfolder '/TempSubjectsMinimaAssignments*'])



