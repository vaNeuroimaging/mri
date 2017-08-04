hem = 'L';
parcelID = [18191];
networksizeminimum = 20;


datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
iscifti = 2;
%hem = 'R';

%kdenval = .2;


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

parcelnum = find(parcelIDs==parcelID);

disp(['Parcel ' num2str(parcelnum)])


Y = pdist(squeeze(parcelcorrelpatterns(:,parcelnum,:)),'correlation');           % 'pdist' converts the square adjacency matrix to a
        %  1 x n matrix so that the function linkage can construct the tree
        
        clustering = linkage(Y, 'average');      % 'linkage' computes the data to construct the tree
        % 'average' refers to the UPGMA algorithm
        
        clusters = cluster(clustering, 'MaxClust', [1:80]);
        
        subject_correlmat = paircorr_mod(squeeze(parcelcorrelpatterns(:,parcelnum,:))');
        
        for numclust = 1:size(clusters,2)
            Qvals(numclust) = M_calc_modularity(clusters(:,numclust),subject_correlmat);
        end
        [maxQval maxQind] = max(Qvals);
        
        [H,T,perm] = dendrogram(clustering, 0, 'orientation','left', 'colorthreshold', clustering(end-maxQind,3));  
        orient landscape; 
        
        cophenetic_r(parcelnum) = cophenet(clustering, Y);


%%



cd(outputfolder)


thisthresh_clusters = clusters(:,maxQind);

communities = unique(thisthresh_clusters);
for comnum = 1:length(communities)
    thisthresh_clusters(thisthresh_clusters==communities(comnum)) = communities(comnum) .* (nnz(thisthresh_clusters==communities(comnum))>=networksizeminimum);
end
communities = unique(thisthresh_clusters);
communities(communities==0) = [];
communityconnections = zeros(length(communities),1);
communitysize = zeros(length(communities),1);
for comnum = 1:length(communities)
    subindices = (thisthresh_clusters==communities(comnum));
    %communitypattern = FisherTransform(squeeze(mean(parcelcorrelpatterns(subindices,parcelnum,:),1)));
    communitypattern = squeeze(mean(parcelcorrelpatterns(subindices,parcelnum,:),1));
    
    output = zeros(32492,1);
    output(logical(mask)) = communitypattern(1:nnz(mask) + + (strcmp(hem,'R') * ncortexLverts));
    save(gifti(single(output)),[outputfolder 'Parcel' num2str(parcelID) '_cluster' num2str(comnum) '_connectivity_' hem '.func.gii'])
    
end




