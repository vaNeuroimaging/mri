hems = {'L','R'};

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    allkdenthreshs = [.1 .15 .2];
    
    for kdenval = allkdenthreshs

datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
iscifti = 2;
%hem = 'R';

%kdenval = .2;
kdenthresholds = [kdenval kdenval];
kdeninterval = .01;
networksizeminimum = 10;

outputfolder = ['/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_allparcels/'];
parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_crossthresh_watershedmerge.func.gii'];
%parcelfilename = ['/data/cn4/evan/ROIs/264_ROIs_' hem '.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;


system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample ' parcelfilename ' /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC ' parcelfilename(1:end-9) '_164.func.gii -largest'])

parcels_upsampled = gifti([parcelfilename(1:end-9) '_164.func.gii']); parcels_upsampled = parcels_upsampled.cdata;

baddata = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); baddata = baddata.cdata;
baddata = baddata<750;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0) =[];

totaltimeseries = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_LR_timeseries.func.gii';
totaltmask = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_total_tmask.txt';

% parcelIDs_incrap = unique(parcels .* baddata);
% for badparcelID = parcelIDs_incrap'
%     parcelIDs(parcelIDs == badparcelID) = [];
% end

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;




if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    parcels = parcels_upsampled;
    ncortexLverts = length(parcels);
end


for parcelnum = 1:length(parcelIDs)
    parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
    
end



[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{1} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{1});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    nsubcortverts = size(subtimecourse,1) - (ncortexLverts+ncortexRverts);


%%

% for subnum = 1:length(subjects)
%     
%     disp(['Subject ' num2str(subnum)])
%     
%     evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' Temp.func.gii']);
%     subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
%     if ~isempty(tmasklist)
%         tmask = load(tmasks{subnum});
%         subtimecourse = subtimecourse(:,logical(tmask));
%     end
%     subtimecourse(isnan(subtimecourse)) = 0;
%     
%     if subnum==1
%         parcelcorrelpatterns = zeros(length(subjects),length(parcelIDs),size(subtimecourse,1));
%         nsubcortverts = size(subtimecourse,1) - (ncortexLverts+ncortexRverts);
%     end
%     
%     parceltimecourses = zeros(length(parcelIDs),size(subtimecourse,2));
%     
%     for parcelnum = 1:length(parcelIDs)
%         parceltimecourses(parcelnum,:) = mean(subtimecourse(parcelindices{parcelnum},:),1);
%         parcelcorrelpatterns(subnum,parcelnum,:) = paircorr_mod(parceltimecourses(parcelnum,:)',subtimecourse');
%         
%     end
%     
% %     parcelcorrelmat = paircorr_mod(parceltimecourses(logical(variableparcelindex),:)',parceltimecourses');
% %     
% %     allsubcorrelmats(:,subnum) = reshape(parcelcorrelmat,numel(parcelcorrelmat),1);
%     
% end
% 
% 
% 
% parcelcorrelpatterns(isnan(parcelcorrelpatterns)) = 0;
%save([outputfolder 'parcelcorrelpatterns_' hem '.mat'],'parcelcorrelpatterns','-v7.3')
load([outputfolder 'parcelcorrelpatterns_' hem '.mat'])
%%
% consensusfile = '/data/cn4/evan/RestingState/Consensus/ConsensusMapvFinal.dtseries.nii';
% evalc(['!wb_command -cifti-convert -to-gifti-ext ' consensusfile ' Temp.func.gii']);
% consensusdata = gifti('Temp.func.gii'); consensusdata = consensusdata.cdata;
% if iscifti==2
%     tempconsensus = zeros(ncortexLverts+ncortexRverts+nsubcortverts,1);
%     
%     giftispace_temp = zeros(32492,1);
%     ciftispace1_maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'); ciftispace1_maskL = ciftispace1_maskL.cdata;
%     giftispace_temp(ciftispace1_maskL==0) = consensusdata(1:nnz(ciftispace1_maskL==0));
%     tempconsensus(1:ncortexLverts) = giftispace_temp(logical(maskL.cdata));
%     
%     giftispace_temp = zeros(32492,1);
%     ciftispace1_maskR = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'); ciftispace1_maskR = ciftispace1_maskR.cdata;
%     giftispace_temp(ciftispace1_maskR==0) = consensusdata(nnz(ciftispace1_maskL==0)+1:nnz(ciftispace1_maskL==0)+nnz(ciftispace1_maskR==0));
%     tempconsensus(ncortexLverts+1:ncortexLverts+ncortexRverts) = giftispace_temp(logical(maskR.cdata));
%     
%     tempconsensus(ncortexLverts+ncortexRverts+1:end) = consensusdata(nnz(ciftispace1_maskL==0)+nnz(ciftispace1_maskR==0)+1:end);
%     
%     consensusdata = tempconsensus;
%     clear tempconsensus
% end
% 
%     
% networkIDs = unique(consensusdata);
% networkIDs(networkIDs==18) = [];
% networkIDs(networkIDs<=0) = [];
% %networkIDs = [1:3,5:17];
% %networkIDs = [1:5,7:16];
% 
% totaltimeseriesdata = gifti(totaltimeseries); totaltimeseriesdata = totaltimeseriesdata.cdata;
% totaltmaskdata = load(totaltmask);
% totaltimeseriesdata = totaltimeseriesdata(:,logical(totaltmaskdata));
% totaltimeseriesdata(isnan(totaltimeseriesdata)) = 0;
% 
% networkcorrelpattern = zeros(size(totaltimeseriesdata,1),length(networkIDs));
% for networknum = 1:length(networkIDs)
%     networktimecourse = mean(totaltimeseriesdata(consensusdata==networkIDs(networknum),:),1);
%     networkcorrelpattern(:,networknum) = paircorr_mod(totaltimeseriesdata',networktimecourse');
% end
% networkcorrelpattern(isnan(networkcorrelpattern)) = 0;

load /data/cn4/evan/RestingState/Consensus/Network_avgconnectivity_xdist20.mat

%surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
%sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
sphere = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii']);
[phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
thetavals = -pi/2 : 1/70 : pi/2;
surf_withlines = zeros(size(phi));
%randvals = randperm(length(thetavals)-1);
reorderedvals = ceil([1:length(thetavals)-1]./2);
addfactor = repmat([0 1],1,ceil(length(reorderedvals)/2)) .* ceil((length(thetavals)-1)/2);
reorderedvals = reorderedvals + addfactor(1:length(reorderedvals));
for i = 1:length(thetavals)-1
    theseindices = intersect(find(theta>thetavals(i)), find(theta>thetavals(i+1)));
    surf_withlines(theseindices) = reorderedvals(i);
end
surf_withlines = surf_withlines+1;
%%
finaloutput = zeros(size(parcels_upsampled));
    
for parcelnum = 1:length(parcelIDs)
    disp(['Parcel ' num2str(parcelnum)])
    
sub_similaritymat = paircorr_mod(squeeze(parcelcorrelpatterns(:,parcelnum,:))');

% imagesc(sub_similaritymat)
% set(gcf,'Nextplot','new');

corrmatname = 'subject_corrmat.mat';
save(corrmatname,'sub_similaritymat');

roifilename = 'subjects.roi';
quickroifile(length(subjects),roifilename);

prmfilename = 'prmfile.txt';
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,['watershed' hem],'-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(1)),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdeninterval),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(2)),'-append','delimiter','');
dlmwrite(prmfilename,outputfolder,'-append','delimiter','');
dlmwrite(prmfilename,num2str(.1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(12),'-append','delimiter','');
dlmwrite(prmfilename,num2str(5),'-append','delimiter','');
dlmwrite(prmfilename,'kden','-append','delimiter','');

%
graphcluster_Evan(prmfilename,'thr',[],1,0,'infomap');

pause(1)

% Modify color assignments with miminum network size criterion

filesinoutputfolder = dir(outputfolder);
for i=3:length(filesinoutputfolder);
    folderdatecreated(i-2) = filesinoutputfolder(i).datenum;
end

[maxval maxindex] = max(folderdatecreated);
trueoutputfolder = [outputfolder '/' filesinoutputfolder(maxindex+2).name];

cd(trueoutputfolder)

simple_assigns = modify_clrfile_nopic('simplify','rawassn.txt',networksizeminimum);

simple_assigns = simple_assigns(:,1);
cd(outputfolder)

%rmdir(trueoutputfolder)


parcel_withlines = surf_withlines .* (parcels_upsampled==parcelIDs(parcelnum));
lines_within_parcel = unique(parcel_withlines);
lines_within_parcel(lines_within_parcel==0)=[];

communities = unique(simple_assigns);
communities(communities==-1) = [];
communityconnections = zeros(length(communities),1);
communitysize = zeros(length(communities),1);
for comnum = 1:length(communities)
    subindices = (simple_assigns==communities(comnum));
    %communitypattern = FisherTransform(squeeze(mean(parcelcorrelpatterns(subindices,parcelnum,:),1)));
    communitypattern = squeeze(mean(parcelcorrelpatterns(subindices,parcelnum,:),1));
    
    correlations_with_network_patterns = paircorr_mod(communitypattern,networkcorrelpattern);
    [ign maxindex] = max(correlations_with_network_patterns);

%     networkcompareTstats = zeros(length(networkIDs),1);
%     for networknum = 1:length(networkIDs)
%         [H,P,CI,STATS] = ttest2(communitypattern(consensusdata==networkIDs(networknum)),communitypattern(consensusdata~=networkIDs(networknum)));
%         networkcompareTstats(networknum) = STATS.tstat;
%     end
%     [ign maxindex] = max(networkcompareTstats);
    communityconnections(comnum) = networkIDs(maxindex);
    communitysize(comnum) = nnz(simple_assigns==communities(comnum));
    
end

cumulativecommunitysize = cumsum(communitysize);

unique_connections = unique(communityconnections);

community_perc = 0;
for comnum = 1:length(unique_connections)
    %community_perc = [community_perc cumulativecommunitysize(comnum)/sum(communitysize)];
    %lineval_indices = [round(community_perc(comnum) .* length(lines_within_parcel))+1 : round(community_perc(comnum+1) .* length(lines_within_parcel))];
    lineval_indices = [(round((comnum-1)/length(unique_connections)*length(lines_within_parcel)) +1) : (round((comnum)/length(unique_connections)*length(lines_within_parcel)))];
    linevals_for_this_community = lines_within_parcel(lineval_indices);
    for lineval = linevals_for_this_community'
        finaloutput(parcel_withlines==lineval) = unique_connections(comnum);
    end
end

end

%save(gifti(single(finaloutput)),['Parcel_networkIDs_subjectvariability_kden_' num2str(kdenthresholds(1)) '_minsize_' num2str(networksizeminimum) '_' hem '_164.func.gii']);
save(gifti(single(finaloutput)),['264_networkIDs_subjectvariability_kden_' num2str(kdenthresholds(1)) '_minsize_' num2str(networksizeminimum) '_' hem '_164.func.gii']);
    end
end