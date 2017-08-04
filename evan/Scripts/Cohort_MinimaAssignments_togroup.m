minimadistances = [2 3 4 5 6 7 8 9 10];
costparameter = 0; %relative weighting of distance vs connectivity pattern
distthreshvals = [10 12 14 16 18 20 22 24 26 28 30 32 34];

hem = 'L';
groupcohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list';
cohortfiles = {'/data/cn4/laumannt/fcMapping_redux/C1_DATALIST_TMASK.txt','/data/cn4/laumannt/fcMapping_redux/C2_DATALIST_TMASK.txt','/data/cn4/laumannt/fcMapping_redux/C3_DATALIST_TMASK.txt'};
cohortedgemaps = {'smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_C1_avg_edge_avg_smooth_L_noalone.func.gii','smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_C2_avg_edge_avg_smooth_L_noalone.func.gii','smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_C3_avg_edge_avg_smooth_L_noalone.func.gii'};
cohortedgemapfolder = '/data/cn4/laumannt/fcMapping_redux/surfsmooth_bandpass_surfreg_32k_gradients/';
edgemapfolder = '/data/cn4/laumannt/left_hem/';
groupedgemapfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/MinimaAssignments/';
groupedgefile = 'smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_AllC_avg_edge_avg_smooth_L_noalone.func.gii';
tmaskname = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
medialmask = '/data/cn4/evan/RestingState/FC_Mapping_120/medialmask.func.gii';
subjBOLDdatafolder = '/data/cn4/laumannt/fcMapping_redux/surfsmooth_bandpass_surfreg_32k_gradients/';

makenewgroupstuff = 0;

outputstatsfile = [outputfolder '/Statistics_cohortmatch_manydist_manydistthresh.txt'];

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
    
    [ign subnum] = system(['cat ' groupcohortfile ' | wc -l']);
    subnum = str2num(subnum);
    
    avgcorrelmat = zeros(length(groupminimavertices),length(groupminimadata));
    
    
    disp('Calculating average connectivity patterns of group minima');
    
    for s = 1:subnum
        
        if makenewgroupstuff
            
            surfBOLDfilename = dir([subjBOLDdatafolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii']);
            surfBOLDfile{s} = [subjBOLDdatafolder '/' surfBOLDfilename(1).name];
            
            tmask{s} = load(tmaskfiles{s});
            
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
    
    
    
    for cohort = 1:length(cohortfiles)
        clear tmask
        
        
        disp(['Cohort ' num2str(cohort)]);
        
        disp('    calculating minima');
        
        cohortminimametric = [cohortedgemaps{cohort}(1:end-9) '_minima_dist' num2str(minimadist) '.func.gii'];
        
        try cohortminimadata = gifti([outputfolder '/' cohortminimametric]);
        catch
            metric_minima([cohortedgemapfolder '/' cohortedgemaps{cohort}],minimadist,[],hem,outputfolder,cohortminimametric);
            cohortminimadata = gifti([outputfolder '/' cohortminimametric]);
        end
        cohortminimadata = cohortminimadata.cdata;
        cohortminimadata(medialmaskind) = 0;
        cohortminimavertices = find(cohortminimadata);
        
        [cohortsubjects cohorttmaskfiles] = textread(cohortfiles{cohort},'%s%s');
        
        disp('    calculating connectivity of minima');
        
        try load([outputfolder '/avgcorrelmat_cohort' num2str(cohort) '_dist' num2str(minimadist)])
        catch
            cohortavgcorrelmat = zeros(length(cohortminimavertices),length(groupminimadata));
            
            for s = 1:length(cohortsubjects)
                
                surfBOLDfilename = dir([subjBOLDdatafolder '/' cohortsubjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii']);
                surfBOLDfile{s} = [subjBOLDdatafolder '/' surfBOLDfilename(1).name];
                
                tmask{s} = load(cohorttmaskfiles{s});
                
                surf_BOLD = load(surfBOLDfile{s});
                surf_BOLD(:,1) = [];
                surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
                
                correlmat = FisherTransform(paircorr_mod(surf_BOLD(cohortminimavertices,:)',surf_BOLD'));
                correlmat(:,medialmaskind) = 0;
                correlmat(logical(isnan(correlmat))) = 0;
                
                cohortavgcorrelmat = cohortavgcorrelmat + correlmat;
                
            end
            
            cohortavgcorrelmat = cohortavgcorrelmat ./ length(cohortsubjects);
            save([outputfolder '/avgcorrelmat_cohort' num2str(cohort) '_dist' num2str(minimadist)],'cohortavgcorrelmat');
        end
        
        disp('    calculating assignments');
        
        cohortvsgroupcorrel = paircorr_mod(avgcorrelmat',cohortavgcorrelmat');
        
        cohortvsgroupdist = distances(:,cohortminimavertices);
        
        for distthreshnum = 1:length(distthreshvals)
            distthresh = distthreshvals(distthreshnum);
            disp(['        distance threshold: ' num2str(distthresh)]);
            
            toofarmatrix = single((cohortvsgroupdist>distthresh));
            toofarmatrix(logical(toofarmatrix)) = Inf;
            
            cohortvsgroupcost = (cohortvsgroupdist ./ maxdist) .* costparameter + toofarmatrix  + (1-(cohortvsgroupcorrel+1)/2);
            
            [assignment, cost] = munkres(cohortvsgroupcost);
%             [assignment, cost, ~, ~, ~] = lapjv(cohortvsgroupcost);
%             if size(cohortvsgroupcost,1)>size(cohortvsgroupcost,2)
%                 newassignment = zeros(size(cohortvsgroupcost,1),1);
%                 for i=1:length(newassignment)
%                     if ~isempty(find(assignment==i));
%                         newassignment(i) = find(assignment==i);
%                     else
%                         newassignment(i) = 0;
%                     end
%                 end
%                 assignment = newassignment;
%             end
            
            for minimanum = 1:length(groupminimavertices)
                if cohort==1; outputmetric(:,minimanum,distthreshnum) = zeros(size(groupminimadata),1); end
                if assignment(minimanum) > 0 && ~isinf(cohortvsgroupcost(minimanum,assignment(minimanum)))
                    outputmetric(cohortminimavertices(assignment(minimanum)),minimanum,distthreshnum) = 1;
                    texttowrite = [num2str(minimadist) '   ' num2str(distthresh) '   Cohort' num2str(cohort) '   ' num2str(groupminimavertices(minimanum)) '   ' num2str(cohortminimavertices(assignment(minimanum))) '   ' num2str(cohortvsgroupcorrel(minimanum,assignment(minimanum))) '   ' num2str(cohortvsgroupdist(minimanum,assignment(minimanum))) '   ' num2str(cohortvsgroupcost(minimanum,assignment(minimanum)))];
                else
                    texttowrite = [num2str(minimadist) '   ' num2str(distthresh) '   Cohort' num2str(cohort) '   ' num2str(groupminimavertices(minimanum)) '   None   None   None   None'];
                end
                outputmetric(groupminimavertices(minimanum),minimanum,distthreshnum) = 30;
                
                dlmwrite(outputstatsfile,texttowrite,'-append','delimiter','');
                
            end
            
        end
        
    end
    
    disp('Writing output');
    
    for distthreshnum = 1:length(distthreshvals)
        distthresh = distthreshvals(distthreshnum);
        
        %save(gifti(single(squeeze(outputmetric(:,:,distthreshnum)))),[outputfolder '/CohortMinimaAssignments_dist' num2str(minimadist) '_param' num2str(costparameter) '_distthresh' num2str(distthresh) '.func.gii'],'ExternalFileBinary');
        
    end
    
end





