minimadistances = [5]; %minima distance parameter of group watershed
overlapcostparameters = [0]; %relative weighting of overlap vs connectivity pattern in the assignment costs
overlapthreshvals = [.00001 .05 .1 .15 .2 .25 .3 .4 .5]; %percent overlap threshold that parcels must exceeed to be considered as a possible match 
connectthresh = []; %minimum similarity of connectivity patterns w/group that parcels must exceeed to be considered as a possible match 
fracmaxh = .7; %watershed erosion parameter (1=not eroded)
subminimadistadjust = 0; %use if you want subjects to have different minima distance parameter than the group's; the subject parameter will be +this factor

hem = 'L';
cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list'; %list of subjects
tmaskname = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt'; %list of tmask files
edgemapfolder = '/data/cn4/laumannt/left_hem_edge/'; %location of subject edge maps
groupedgemapfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/WatershedAssignments/'; %location of group edge map
groupedgefile = 'fUC_smooth_smooth2.55.func.gii'; %name of group edge map
subjBOLDdatafolder ='/data/cn4/laumannt/fcMapping_redux/surfsmooth_bandpass_surfreg_32k_gradients/'; % location of subject surface timecourses '/data/cn4/evan/RestingState/FC_Mapping_120/WatershedAssignments/';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/WatershedAssignments/'; %location watersheds and match data will be written to
medialmask = '/data/cn4/evan/RestingState/FC_Mapping_120/medialmask.func.gii'; %medial mask file
gooddatafile = '/data/cn4/evan/fsaverage_LR32k/all_meanimage_L.func.gii'; %file indicating signal dropout; only parcels not within these regions will be considered for a match

surfacecoordfile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.coord.gii';
surfacetopofile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.32k_fs_LR.topo.gii';

outputstatsfile = [outputfolder '/WatershedAssignments_minimadist5_manyoverlapthresh_comparegroupparcel_localexclusion.txt']; %match quality data writes to this file

delete(outputstatsfile);

%open the output stats file for writing and write a header
fid = fopen(outputstatsfile,'at'); 
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\n\r\','Minimadist','Subminimadist','Overlapthresh','OverlapCost','Subject','GroupWatershed','SubjectWatershed','Correl','Overlap'); %write the output file header
fclose(fid);
dlmwrite(outputstatsfile,' ','-append');

%load medial mask
medialmaskdata = gifti(medialmask);
medialmaskind = find(~medialmaskdata.cdata);
notmedialmaskind = find(medialmaskdata.cdata);

thisdir = pwd;

%loop through minima distance parameters
for minnum = 1:length(minimadistances)
    minimadist = minimadistances(minnum);
    subminimadist = minimadist + subminimadistadjust;

    groupminimafile = [groupedgefile(1:end-9) '_minima_dist' num2str(minimadist) '.func.gii'];
    groupwatershedfile = [groupminimafile(1:end-9) 'watershed_rois.func.gii'];
    
    %try to load group watersheds. If there aren't any, make them.
    try groupwatershed = gifti_evan([outputfolder '/' groupwatershedfile]); disp('Group watersheds loaded')
    catch
        disp('Calculating group watersheds')
        metric_minima([groupedgemapfolder '/' groupedgefile],minimadist,[],hem,outputfolder,groupminimafile);
        watershed_algorithm_Evan([groupedgemapfolder '/' groupedgefile],[outputfolder '/' groupminimafile],200,fracmaxh,outputfolder,groupminimafile(1:end-9));
        groupwatershed = gifti_evan([outputfolder '/' groupwatershedfile]);
    end
    
groupwatershed = groupwatershed.data{1}.data;

%load the file indicating signal dropout and figure out where the dropout is
gooddata = gifti(gooddatafile);
gooddata = gooddata.cdata;
baddatamask = gooddata < 750;

%figure out which group parcels we care about (that aren't in signal dropout areas) 
groupwatershed_inbaddata = groupwatershed .* repmat(baddatamask,[1 size(groupwatershed,2)]);
goodparcels = find(sum(groupwatershed_inbaddata,1)<=5);

%restrict our watersheds to those parcels
groupwatershed = groupwatershed(:,goodparcels);

%load the list of subjects
[subjects tmaskfiles]=textread(tmaskname,'%s%s');
subnum = length(subjects);

%try to load the average connectivity patterns of the group watersheds. If these don't exist, prep the variables needed to make them 
makenewgroupstuff = 0;
try load([outputfolder '/avgcorrelmat_' num2str(minimadist)]); disp('Average connectivity patterns of group watersheds loaded');
catch
    makenewgroupstuff = 1;
    disp('Calculating average connectivity patterns of group watersheds');
    avgcorrelmat = zeros(size(groupwatershed,2),size(groupwatershed,1));
end

%subject loop
for s = 1:subnum
    
     
    surfBOLDfilename = [subjBOLDdatafolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii'];
    %check whether the "noHEAD" version of this subject's surface tiomecourse exists. If not, make it
    if ~exist(surfBOLDfilename)
        if ~exist([outputfolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii'])
            evalc(['!awk ''NF > 25'' ' surfBOLDfilename(1:end-16) '.func.gii > ' outputfolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii']);
        end
        surfBOLDfilename = [outputfolder '/' subjects{s} '_time_smooth2.55_L_32k_fsLR_noHEAD.func.gii'];
    end
    surfBOLDfile{s} = surfBOLDfilename;
    
    
    tmask{s} = load(tmaskfiles{s});
    
    %load this subject's data and add it to the average connectivity patterns of the group watersheds
    if makenewgroupstuff
        disp(['    Subject #' num2str(s) ': ' subjects{s}])
        surf_BOLD = load(surfBOLDfile{s});
            
        surf_BOLD(:,1) = [];
        surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
        
        %loop through group watershed parcels and get this subject's connectivity pattern for each 
        for watershednum = 1:size(groupwatershed,2)
        
            correlmat(watershednum,:) = FisherTransform(paircorr_mod(mean(surf_BOLD(logical(groupwatershed(:,watershednum)),:),1)',surf_BOLD'));
        
        end
        
        correlmat(:,medialmaskind) = 0;
        correlmat(logical(isnan(correlmat))) = 0;
        
        %add this subject's connectivity pattern to the average 
        avgcorrelmat = avgcorrelmat + correlmat;
    end
    
end

%finalize and save the average connectivity patterns
if makenewgroupstuff
    avgcorrelmat = avgcorrelmat ./ subnum;
    save([outputfolder '/avgcorrelmat_' num2str(minimadist)],'avgcorrelmat');
end

%set up variables to check the connectivity of the parcel itself 
indmatchvalues_ingroupparcels = zeros(size(avgcorrelmat,2),subnum);
groupmatchvalues = zeros(size(avgcorrelmat,2),subnum);


disp('Conducting watershed assignments');
for s = 1:subnum
    
    disp(['Subject #' num2str(s) ': ' subjects{s}])
    subjectedgemetric = dir([edgemapfolder '/' subjects{s} '_avg_edge_avg_smooth_' hem '_noalone_smooth2.55.func.gii']);
    subjectedgemetric = subjectedgemetric(1).name;
    
    subjectminimametric = [subjectedgemetric(1:end-9) '_minimadist' num2str(subminimadist) '.func.gii'];
    subjectwatershedmetric = [subjectminimametric(1:end-9) 'watershed_rois.func.gii'];
    
    %try to load this subject's watershed. If that doesn't work, make it.
    try subjectwatershed = gifti_evan([outputfolder '/' subjectwatershedmetric]); disp('    watersheds loaded')
    catch
        disp('    generating watersheds')
        fprintf('     ')
        metric_minima([edgemapfolder '/' subjectedgemetric],subminimadist,[],hem,outputfolder,subjectminimametric);
        watershed_algorithm_Evan([edgemapfolder '/' subjectedgemetric],[outputfolder '/' subjectminimametric],200,fracmaxh,outputfolder,subjectminimametric(1:end-9));
        subjectwatershed = gifti_evan([outputfolder '/' subjectwatershedmetric]);
    end
    
    subjectwatershed = subjectwatershed.data{1}.data;
    
    %delete parcels in the medial mask
    subjectwatershed(medialmaskind,:) = 0;
    
    disp('    calculating connectivity of watersheds')
    
    %load this subject's surface timecourse
    surf_BOLD = load(surfBOLDfile{s});
    surf_BOLD(:,1) = [];
    surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
    
    %for each of this subject's watershed parcels, calculate the connectivity pattern and the overlap with all group parcels 
    for watershednum = 1:size(subjectwatershed,2)
        
        subjcorrelmat(watershednum,:) = FisherTransform(paircorr_mod(mean(surf_BOLD(logical(subjectwatershed(:,watershednum)),:),1)',surf_BOLD'));
        
        overlap = repmat(subjectwatershed(:,watershednum),[1 size(groupwatershed,2)]) .* groupwatershed;
        
        percentoverlap(:,watershednum) = mean([(sum(overlap,1) ./ sum(groupwatershed,1)); (sum(overlap,1) ./ sum(subjectwatershed(:,watershednum)))],1)';
        
    end
    
    %remove NaNs and values in the medial mask
    subjcorrelmat(:,medialmaskind) = 0;
    subjcorrelmat(logical(isnan(subjcorrelmat))) = 0;
    
    
    %figure out subject connectivity patterns for group parcels
    for parcelnum = 1:size(groupwatershed,2)
        
            string{parcelnum} = ['    calculating subject connectivity of group parcel ' num2str(parcelnum) ' of ' num2str(size(groupwatershed,2))];
            if parcelnum==1; fprintf('%s',string{parcelnum}); else fprintf([repmat('\b',1,length(string{parcelnum-1})) '%s'],string{parcelnum}); end
        
            %figure out this subject's connectivity pattern for this group parcel
            watershedmap = FisherTransform(paircorr_mod(mean(surf_BOLD(logical(groupwatershed(:,parcelnum)),:),1)',surf_BOLD'));
            watershedmap(logical(isnan(watershedmap))) = 0;
            
            %figure out which nodes we want to exclude from the comparison; i.e., the nodes actually inside the group parcel
            nodestocorrelate = find(groupwatershed(:,parcelnum)==0 & medialmaskdata.cdata);
            
            %correlated the subject correlation map against the group correlation map
            correlation = paircorr_mod(watershedmap(nodestocorrelate)',avgcorrelmat(parcelnum,nodestocorrelate)');
        
            %write the resulting correlation value
            texttowrite = [num2str(minimadist) '   Group   Group   Group   ' subjects{s} '   ' num2str(goodparcels(parcelnum)) '   Group   ' num2str(correlation) '   1'];
            dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
            
            %save the resulting correlation value to be written to an image 
            groupmatchvalues(logical(groupwatershed(:,parcelnum)),s) = correlation;
            
    end
    
    disp(' ')
    
    %figure out which watershed pairs have at least a little overlap
    atleastsomeoverlapmatrix = single(percentoverlap > 0);
    [groupoverlapindices, subjoverlapindices] = find(atleastsomeoverlapmatrix);
    
    subvsgroupcorrel = ones(size(groupwatershed,2),size(subjectwatershed,2)) * Inf;
    
     
    %for each pair of group vs subject parcels that overlap at least a little
    for parcelcombonum = 1:length(groupoverlapindices) %1:size(groupwatershed,2)
        
                
                string{parcelcombonum} = ['    calculating subject-to-group parcel match costs for comparison ' num2str(parcelcombonum) ' of ' num2str(length(groupoverlapindices))];
                if parcelcombonum==1; fprintf('%s',string{parcelcombonum}); else fprintf([repmat('\b',1,length(string{parcelcombonum-1})) '%s'],string{parcelcombonum}); end
                
                %figure out which nodes we want in our comparison
                nodestocorrelate = find(subjectwatershed(:,subjoverlapindices(parcelcombonum))==0 & groupwatershed(:,groupoverlapindices(parcelcombonum))==0 & medialmaskdata.cdata);
                
                %compare the group parcel's group correlation map vs the subject parcel's correlation map in that subject
                subvsgroupcorrel(groupoverlapindices(parcelcombonum),subjoverlapindices(parcelcombonum)) = paircorr_mod(single(avgcorrelmat(groupoverlapindices(parcelcombonum),nodestocorrelate))',single(subjcorrelmat(subjoverlapindices(parcelcombonum),nodestocorrelate))');
                
                
     end
        
     
    disp(' ')
    disp('    calculating assignments');
    
    %for each overlap threshold we're testing
    for overlapthreshnum = 1:length(overlapthreshvals)
        disp(['        overlap threshold: ' num2str(overlapthreshvals(overlapthreshnum))]);
        
        for overlapcostnum = 1:length(overlapcostparameters)
            
            %disp(['          overlap cost parameter: ' num2str(overlapcostparameters(overlapcostnum))]);
        
        %figure out which parcel pairs overlap too little
        toolittleoverlapmatrix = single(percentoverlap < overlapthreshvals(overlapthreshnum));
        %and set those pairs to infinity (to be added to the cost function below) 
        toolittleoverlapmatrix(logical(toolittleoverlapmatrix)) = Inf;
        
        %calculate the cost function for all parcel pairs 
        subjvsgroupcost =  (1-(subvsgroupcorrel+1)/2) + (1-percentoverlap)*overlapcostparameters(overlapcostnum) + toolittleoverlapmatrix;
        
        %calculate the assignments
        [assignment, cost] = munkres(subjvsgroupcost);

        %for each group parcel
        for parcelnum = 1:size(groupwatershed,2)
            
            %set up an output image
            if s==1; outputmetric(:,parcelnum,overlapthreshnum,overlapcostnum) = zeros(size(groupwatershed,1),1); end
            
            %if this group parcel has any assignment
            if assignment(parcelnum) > 0 && ~isinf(subjvsgroupcost(parcelnum,assignment(parcelnum)))
                
                %save the match quality in one output image
                indmatchvalues_ingroupparcels(logical(groupwatershed(:,parcelnum)),s,overlapthreshnum) = subvsgroupcorrel(parcelnum,assignment(parcelnum));
                
                %save the subject parcel location in another output image
                outputmetric(:,parcelnum,overlapthreshnum,overlapcostnum) = outputmetric(:,parcelnum,overlapthreshnum,overlapcostnum) + subjectwatershed(:,assignment(parcelnum));
                
                %set the statistics to write out
                texttowrite = [num2str(minimadist) '   ' num2str(subminimadist) '   ' num2str(overlapthreshvals(overlapthreshnum)) '   ' num2str(overlapcostparameters(overlapcostnum)) '   ' subjects{s} '   ' num2str(goodparcels(parcelnum)) '   ' num2str(assignment(parcelnum)) '   ' num2str(subvsgroupcorrel(parcelnum,assignment(parcelnum))) '   ' num2str(percentoverlap(parcelnum,assignment(parcelnum)))];
            else
                %set up to write out the lack of a match
                texttowrite = [num2str(minimadist) '   ' num2str(subminimadist) '   ' num2str(overlapthreshvals(overlapthreshnum)) '   ' num2str(overlapcostparameters(overlapcostnum)) '   ' subjects{s} '   ' num2str(goodparcels(parcelnum))];
            end
            
           %write out the statistics
           dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
            
        end
        end
    end
    
    clear subjcorrelmat percentoverlap subvsgroupcorrel
    
end

disp('Writing output');

groupmatchvalues(logical(isnan(groupmatchvalues))) = 0;

%save the average and std dev of the correlations only within the group parcels 
save(gifti(single(squeeze(mean(groupmatchvalues,2)))),[outputfolder '/AvgGroupWatershedCorrel_dist' num2str(minimadist) '.func.gii']);
save(gifti(single(squeeze(std(groupmatchvalues,[],2)))),[outputfolder '/StdGroupWatershedCorrel_dist' num2str(minimadist) '.func.gii'])

for overlapthreshnum = 1:length(overlapthreshvals)
    overlapthresh = overlapthreshvals(overlapthreshnum);
    
    %save the average and standard deviation of the sub vs group parcel match quality 
    save(gifti(single(squeeze(mean(indmatchvalues_ingroupparcels(:,:,overlapthreshnum),2)))),[outputfolder '/AvgSubtoGroupMatchWatershedCorrel_subdist' num2str(subminimadist) '_overlapthresh' num2str(overlapthresh) '.func.gii']);
    save(gifti(single(squeeze(std(indmatchvalues_ingroupparcels(:,:,overlapthreshnum),[],2)))),[outputfolder '/StdSubtoGroupMatchWatershedCorrel_subdist' num2str(subminimadist) '_overlapthresh' num2str(overlapthresh) '.func.gii']);

    for overlapcostnum = 1:length(overlapcostparameters)
        overlapcost = overlapcostparameters(overlapcostnum);
        
        %save the location of the subject match
        save(gifti(single(squeeze(outputmetric(:,:,overlapthreshnum,overlapcostnum)))),[outputfolder '/SubjectWatershedAssignments_subdist' num2str(subminimadist) '_overlapthresh' num2str(overlapthresh) '_overlapcost' num2str(overlapcost) '.func.gii'],'ExternalFileBinary');
    
    end
end
clear outputmetric

cd(thisdir)

end





