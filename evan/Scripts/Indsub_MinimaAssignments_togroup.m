minimadistances = [2];
costparameter = 0; %relative weighting of distance vs connectivity pattern
distthreshvals = [10 12 14 16 18 20 22 24 26 28 30 32 34];
connectthresh = [];

hem = 'L';
cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list';
edgemapfolder = '/data/cn4/laumannt/left_hem_edge/';
groupedgemapfolder = '/data/cn4/laumannt/assign_indiv_group/';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/MinimaAssignments/';
groupedgefile = 'AllC_edge_smooth2.55.func.gii';
tmaskname = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
medialmask = '/data/cn4/evan/RestingState/FC_Mapping_120/medialmask.func.gii';
subjBOLDdatafolder = '/data/cn4/laumannt/fcMapping_redux/surfsmooth_bandpass_surfreg_32k_gradients/';

makenewgroupstuff = 0;

outputstatsfile = [outputfolder '/Statistics_manydist_manydistthresh_jv2.txt'];

delete(outputstatsfile);
fid = fopen(outputstatsfile,'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\n\r\','Minimadist','Distthresh','Subject','GroupNode','SubjectNode','Correl','Distance','OverallCost'); %write the output file header
fclose(fid);
dlmwrite(outputstatsfile,' ','-append');

surfacecoordfile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.coord.gii';
surfacetopofile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.32k_fs_LR.topo.gii';

medialmaskdata = gifti(medialmask);
medialmaskind = find(~medialmaskdata.cdata);


for minnum = 1:length(minimadistances)
    minimadist = minimadistances(minnum);

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
    save([outputfolder '/distances_' num2str(minimadist)],'distances')
    system('rm temp.func.gii');
else
    load([outputfolder '/distances_' num2str(minimadist)]);
end

maxdist = max(max(distances));

[subjects tmaskfiles]=textread(tmaskname,'%s%s');

[ign subnum] = system(['cat ' cohortfile ' | wc -l']);
subnum = str2num(subnum);

avgcorrelmat = zeros(length(groupminimavertices),length(groupminimadata));


disp('Calculating average connectivity patterns of group minima');

for s = 1:subnum
%     [ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
%     [ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
%     [ign funcvols{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1}''']);
%     [ign surfdirs{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' cohortfile]);
    
    surfBOLDfilename = dir([subjBOLDdatafolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii']);
    surfBOLDfile{s} = [subjBOLDdatafolder '/' surfBOLDfilename(1).name];
    
    tmask{s} = load(tmaskfiles{s});
    
    if makenewgroupstuff
        disp(['    Subject ' subjects{s}])
        
        surf_BOLD = load(surfBOLDfile{s});
        surf_BOLD(:,1) = [];
        surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
        
        correlmat = FisherTransform(paircorr_mod(surf_BOLD(groupminimavertices,:)',surf_BOLD'));
        correlmat(:,medialmaskind) = 0;
        correlmat(logical(isnan(correlmat))) = 0;
        
        avgcorrelmat = avgcorrelmat + correlmat;
    end
    
end

if makenewgroupstuff
    avgcorrelmat = avgcorrelmat ./ subnum;
    save([outputfolder '/avgcorrelmat_' num2str(minimadist)],'avgcorrelmat');
else
    load([outputfolder '/avgcorrelmat_' num2str(minimadist)]);
end




disp('Conducting minima assignments');
for s = 1:subnum
    
    disp(['Subject ' subjects{s} ':']);
    disp('    calculating minima');
    subjectedgemetric = dir([edgemapfolder '/' subjects{s} '_avg_edge_avg_smooth_' hem '_noalone_smooth2.55.func.gii']);
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
    
    surf_BOLD = load(surfBOLDfile{s});
    surf_BOLD(:,1) = [];
    surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
    
    subjcorrelmat = FisherTransform(paircorr_mod(surf_BOLD(subjectminimavertices,:)',surf_BOLD'));
    subjcorrelmat(:,medialmaskind) = 0;
    subjcorrelmat(logical(isnan(subjcorrelmat))) = 0;
    
    
    disp('    calculating assignments');
    
    subvsgroupcorrel = paircorr_mod(avgcorrelmat',subjcorrelmat');
    
    subjvsgroupdist = distances(:,subjectminimavertices);
    
    if ~isempty(connectthresh)
        toolowmatrix = single(subvsgroupcorrel<connectthresh);
        toolowmatrix(logical(toolowmatrix)) = Inf;
    else
        toolowmatrix = zeros(size(subvsgroupcorrel));
    end
    
    for distthreshnum = 1:length(distthreshvals)
            distthresh = distthreshvals(distthreshnum);
            disp(['        distance threshold: ' num2str(distthresh)]);
            
    toofarmatrix = single((subjvsgroupdist>distthresh));
    toofarmatrix(logical(toofarmatrix)) = Inf;
    
    subjvsgroupcost = (subjvsgroupdist ./ maxdist) .* costparameter + toofarmatrix + toolowmatrix + (1-(subvsgroupcorrel+1)/2);
    
%    [assignment, cost] = munkres(subjvsgroupcost);
    [assignment, cost, ~, ~, ~] = lapjv(subjvsgroupcost);
    if size(subjvsgroupcost,1)>size(subjvsgroupcost,2)
        newassignment = zeros(size(subjvsgroupcost,1),1);
        for i=1:length(newassignment)
            if ~isempty(find(assignment==i));
                newassignment(i) = find(assignment==i);
            else
                newassignment(i) = 0;
            end
        end
        assignment = newassignment;
    end
            
            
    for minimanum = 1:length(groupminimavertices)
        if s==1; outputmetric(:,minimanum,distthreshnum) = zeros(size(groupminimadata),1); end
        if assignment(minimanum) > 0 && ~isinf(subjvsgroupcost(minimanum,assignment(minimanum)))
        outputmetric(subjectminimavertices(assignment(minimanum)),minimanum,distthreshnum) = 1;
        texttowrite = [num2str(minimadist) '   ' num2str(distthresh) '   ' subjects{s} '   ' num2str(groupminimavertices(minimanum)) '   ' num2str(subjectminimavertices(assignment(minimanum))) '   ' num2str(subvsgroupcorrel(minimanum,assignment(minimanum))) '   ' num2str(subjvsgroupdist(minimanum,assignment(minimanum))) '   ' num2str(subjvsgroupcost(minimanum,assignment(minimanum)))];
        else
            texttowrite = [num2str(minimadist) '   ' num2str(distthresh) '   ' subjects{s} '   ' num2str(groupminimavertices(minimanum)) '   None   None   None   None'];
        end 
        outputmetric(groupminimavertices(minimanum),minimanum,distthreshnum) = 30;
        
        dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
                
    end
    end
    
end

disp('Writing output');

for distthreshnum = 1:length(distthreshvals)
        distthresh = distthreshvals(distthreshnum);
        
save(gifti(single(squeeze(outputmetric(:,:,distthreshnum)))),[outputfolder '/SubjectsMinimaAssignments_dist' num2str(minimadist) '_param' num2str(costparameter) '_distthresh' num2str(distthresh) '_lapjv.func.gii'],'ExternalFileBinary');
end
clear outputmetric


end





