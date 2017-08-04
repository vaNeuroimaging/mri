minimadistances = [5];
overlapcostparameters = [0 .2 .4 .6 .8 1]; %relative weighting of overlap vs connectivity pattern
overlapthreshvals = [.2];% .00001 .05 .1 .15 .2 .25 .3 .4 .5];
connectthresh = [];

hem = 'L';
cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list';
edgemapfolder = '/data/cn4/laumannt/left_hem_edge/';
groupedgemapfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/WatershedAssignments/';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/WatershedAssignments/';
groupedgefile = 'fUC_smooth_smooth2.55.func.gii';
tmaskname = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
medialmask = '/data/cn4/evan/RestingState/FC_Mapping_120/medialmask.func.gii';
subjBOLDdatafolder = '/data/cn4/laumannt/fcMapping_redux/surfsmooth_bandpass_surfreg_32k_gradients/';

%makenewgroupstuff = 1;

outputstatsfile = [outputfolder '/WatershedAssignments_minimadist5_overlapthresh.2_manyoverlapcost.txt'];

delete(outputstatsfile);
fid = fopen(outputstatsfile,'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\n\r\','Minimadist','Overlapthresh','OverlapCost','Subject','GroupWatershed','SubjectWatershed','Correl','Overlap'); %write the output file header
fclose(fid);
dlmwrite(outputstatsfile,' ','-append');

surfacecoordfile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.coord.gii';
surfacetopofile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.32k_fs_LR.topo.gii';

medialmaskdata = gifti(medialmask);
medialmaskind = find(~medialmaskdata.cdata);

thisdir = pwd;

for minnum = 1:length(minimadistances)
    minimadist = minimadistances(minnum);

    groupminimafile = [groupedgefile(1:end-9) '_minima_dist' num2str(minimadist) '.func.gii'];
    groupwatershedfile = [groupminimafile(1:end-9) 'watershed_rois.func.gii'];

    try groupwatershed = gifti_evan(groupwatershedfile);
    catch
        
        disp('Calculating group watersheds')
        metric_minima([groupedgemapfolder '/' groupedgefile],minimadist,[],hem,outputfolder,groupminimafile);
        watershed_algorithm_Evan([groupedgemapfolder '/' groupedgefile],[outputfolder '/' groupminimafile],200,1,outputfolder,groupminimafile(1:end-9));
        groupwatershed = gifti_evan(groupwatershedfile);
    
    end
    
groupwatershed = groupwatershed.data{1}.data;


[subjects tmaskfiles]=textread(tmaskname,'%s%s');

%[ign subnum] = system(['cat ' cohortfile ' | wc -l']);
subnum = length(subjects);

avgcorrelmat = zeros(size(groupwatershed,2),size(groupwatershed,1));


disp('Calculating average connectivity patterns of group minima');

makenewgroupstuff = 0;
try load([outputfolder '/avgcorrelmat_' num2str(minimadist)]);
catch
    makenewgroupstuff = 1;
end

for s = 1:subnum
    
    surfBOLDfilename = dir([subjBOLDdatafolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii']);
    %surfBOLDfilename = dir([subjBOLDdatafolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_tmasked.func.gii']);
    surfBOLDfile{s} = [subjBOLDdatafolder '/' surfBOLDfilename(1).name];
    
    tmask{s} = load(tmaskfiles{s});
    
    if makenewgroupstuff
        disp(['    Subject #' num2str(s) ': ' subjects{s}])
        try surf_BOLD = load(surfBOLDfile{s});
        catch
            copyfile([surfBOLDfile{s}(1:end-16) '.func.gii'],[outputfolder '/' surfBOLDfilename(1).name(1:end-16) '.func.gii']);
            altsurfBOLDfile{s} = [outputfolder '/' surfBOLDfilename(1).name];
            system(['awk ''NF > 25'' ' outputfolder '/' surfBOLDfilename(1).name(1:end-16) '.func.gii > ' altsurfBOLDfile{s}])
            surf_BOLD = load(altsurfBOLDfile{s});
        end
            
        surf_BOLD(:,1) = [];
        surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
        
        for watershednum = 1:size(groupwatershed,2)
        
            correlmat(watershednum,:) = FisherTransform(paircorr_mod(mean(surf_BOLD(logical(groupwatershed(:,watershednum)),:),1)',surf_BOLD'));
        
        end
        
        correlmat(:,medialmaskind) = 0;
        correlmat(logical(isnan(correlmat))) = 0;
        
        avgcorrelmat = avgcorrelmat + correlmat;
    end
    
end

if makenewgroupstuff
    avgcorrelmat = avgcorrelmat ./ subnum;
    save([outputfolder '/avgcorrelmat_' num2str(minimadist)],'avgcorrelmat');
end




disp('Conducting watershed assignments');
for s = 1:subnum
    
    disp(['Subject #' num2str(s) ': ' subjects{s}])
    disp('    generating watersheds');
    subjectedgemetric = dir([edgemapfolder '/' subjects{s} '_avg_edge_avg_smooth_' hem '_noalone_smooth2.55.func.gii']);
    subjectedgemetric = subjectedgemetric(1).name;
    
    subjectminimametric = [subjectedgemetric(1:end-9) '_minimadist' num2str(minimadist) '.func.gii'];
    subjectwatershedmetric = [subjectminimametric(1:end-9) 'watershed_rois.func.gii'];
    
    try subjectwatershed = gifti_evan([outputfolder '/' subjectwatershedmetric]);
    catch
        %         temp = gifti([edgemapfolder '/' subjectedgemetric]);
        %         [ign tempmetric] = upper_completion(temp.cdata);
        %         save(gifti(tempmetric),[outputfolder '/' subjectedgemetric(1:end-9) '_uc.func.gii']);
        fprintf('     ')
        metric_minima([edgemapfolder '/' subjectedgemetric],minimadist,[],hem,outputfolder,subjectminimametric);
        watershed_algorithm_Evan([edgemapfolder '/' subjectedgemetric],[outputfolder '/' subjectminimametric],200,1,outputfolder,subjectminimametric(1:end-9));
        subjectwatershed = gifti_evan([outputfolder '/' subjectwatershedmetric]);
    end
    
    subjectwatershed = subjectwatershed.data{1}.data;
    
    %subjectminimadata = subjectminimadata.cdata;
    subjectwatershed(medialmaskind,:) = 0;
    %subjectminimavertices = find(subjectminimadata);
    
    disp('    calculating connectivity of watersheds');
    
    try surf_BOLD = load(surfBOLDfile{s});
    catch
        surf_BOLD = load(altsurfBOLDfile{s});
    end
    surf_BOLD(:,1) = [];
    surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
    
    for watershednum = 1:size(subjectwatershed,2)
        
        subjcorrelmat(watershednum,:) = FisherTransform(paircorr_mod(mean(surf_BOLD(logical(subjectwatershed(:,watershednum)),:),1)',surf_BOLD'));
        
        overlap = repmat(subjectwatershed(:,watershednum),[1 size(groupwatershed,2)]) .* groupwatershed;
        
        percentoverlap(:,watershednum) = mean([(sum(overlap,1) ./ sum(groupwatershed,1)); (sum(overlap,1) ./ sum(subjectwatershed(:,watershednum)))],1)';
        
    end
    
    %subjcorrelmat = FisherTransform(paircorr_mod(surf_BOLD(subjectminimavertices,:)',surf_BOLD'));
    subjcorrelmat(:,medialmaskind) = 0;
    subjcorrelmat(logical(isnan(subjcorrelmat))) = 0;
    
    
    disp('    calculating assignments');
    
    subvsgroupcorrel = paircorr_mod(avgcorrelmat',subjcorrelmat');
    
    %     subjvsgroupdist = distances(:,subjectminimavertices);
    %
    %     if ~isempty(connectthresh)
    %         toolowmatrix = single(subvsgroupcorrel<connectthresh);
    %         toolowmatrix(logical(toolowmatrix)) = Inf;
    %     else
    %         toolowmatrix = zeros(size(subvsgroupcorrel));
    %     end
    %
    for overlapthreshnum = 1:length(overlapthreshvals)
        %             distthresh = distthreshvals(distthreshnum);
        disp(['        overlap threshold: ' num2str(overlapthreshvals(overlapthreshnum))]);
        
        for overlapcostnum = 1:length(overlapcostparameters)
            
            disp(['          overlap cost parameter: ' num2str(overlapcostparameters(overlapcostnum))]);
        
        toolittleoverlapmatrix = single(percentoverlap < overlapthreshvals(overlapthreshnum));
        toolittleoverlapmatrix(logical(toolittleoverlapmatrix)) = Inf;
        
        subjvsgroupcost =  (1-(subvsgroupcorrel+1)/2) + (1-percentoverlap)*overlapcostparameters(overlapcostnum) + toolittleoverlapmatrix;
        
        [assignment, cost] = munkres(subjvsgroupcost);
        %    [assignment, cost, ~, ~, ~] = lapjv(subjvsgroupcost);
        %     if size(subjvsgroupcost,1)>size(subjvsgroupcost,2)
        %         newassignment = zeros(size(subjvsgroupcost,1),1);
        %         for i=1:length(newassignment)
        %             if ~isempty(find(assignment==i));
        %                 newassignment(i) = find(assignment==i);
        %             else
        %                 newassignment(i) = 0;
        %             end
        %         end
        %         assignment = newassignment;
        %     end
        
        
        for parcelnum = 1:size(groupwatershed,2)
            if s==1; outputmetric(:,parcelnum,overlapthreshnum,overlapcostnum) = zeros(size(groupwatershed),1); end
            
            if assignment(parcelnum) > 0 && ~isinf(subjvsgroupcost(parcelnum,assignment(parcelnum)))
                
                outputmetric(:,parcelnum,overlapthreshnum,overlapcostnum) = outputmetric(:,parcelnum,overlapthreshnum,overlapcostnum) + subjectwatershed(:,assignment(parcelnum));
                texttowrite = [num2str(minimadist) '   ' num2str(overlapthreshvals(overlapthreshnum)) '   ' num2str(overlapcostparameters(overlapcostnum)) '   ' subjects{s} '   ' num2str(parcelnum) '   ' num2str(assignment(parcelnum)) '   ' num2str(subvsgroupcorrel(parcelnum,assignment(parcelnum))) '   ' num2str(percentoverlap(parcelnum,assignment(parcelnum)))];
            else
                texttowrite = [num2str(minimadist) '   ' num2str(overlapthreshvals(overlapthreshnum)) '   ' num2str(overlapcostparameters(overlapcostnum)) '   ' subjects{s} '   ' num2str(parcelnum) '   None   None   None'];
            end
            %outputmetric(groupminimavertices(minimanum),minimanum,distthreshnum) = 30;
            
            dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
            
        end
        end
    end
    
    clear subjcorrelmat percentoverlap
    
end

disp('Writing output');

for overlapthreshnum = 1:length(overlapthreshvals)
    overlapthresh = overlapthreshvals(overlapthreshnum);
    for overlapcostnum = 1:length(overlapcostparameters)
        overlapcost = overlapcostparameters(overlapcostnum);
        
        save(gifti(single(squeeze(outputmetric(:,:,overlapthreshnum,overlapcostnum)))),[outputfolder '/SubjectsMinimaAssignments_dist' num2str(minimadist) '_overlapthresh' num2str(overlapthresh) '._overlapcost' num2str(overlapcost) 'func.gii'],'ExternalFileBinary');
    
    end
end
clear outputmetric

cd(thisdir)

end





