%% Cluster Movement Metrics

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');

mask{1} = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'); mask{1} = ~mask{1}.cdata;
mask{2} = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'); mask{2} = ~mask{2}.cdata;

% verts{1} = gifti('/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii'); verts{1} = verts{1}.vertices;
% verts{2} = gifti('/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii'); verts{2} = verts{2}.vertices;


nsurfverts = nnz(mask{1}) + nnz(mask{2});

neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');

pct_overlap = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
size_change = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
centroid_movement = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
border_movement = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
border_movement_directed = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;

load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
distances = distances(1:nsurfverts,1:nsurfverts);

for s = 1:size(networkconnection_bysub,2)
    disp(['Subject ' num2str(s)])
    fprintf('   ')
    [groupmatches, submatches] = networkClusters_matchtogroup(networkconnection_bysub(:,s),distances,mostcommon);
    
    SA_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.midthickness.32k_fs_LR_surfaceareas.func.gii']);
    SA_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.midthickness.32k_fs_LR_surfaceareas.func.gii']);
    SA_LR = [SA_L.cdata(logical(mask{1})) ; SA_R.cdata(logical(mask{2}))];
    
    for i = 1:size(groupmatches,2)
        
        group_borderinds = [];
        groupinds = find(groupmatches(:,i));
        for ind = groupinds(:)'
            indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
            if any(groupmatches(indneighs,i)==0)
                group_borderinds(end+1) = ind;
            end
        end
        
        sub_borderinds = [];
        subinds = find(submatches(:,i));
        for ind = subinds(:)'
            indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
            if any(submatches(indneighs,i)==0)
                sub_borderinds(end+1) = ind;
            end
        end
        
        [mindist, mini] = min(distances(group_borderinds,sub_borderinds),[],2);
        border_movement(group_borderinds,s) = mindist;
        
        mindist_directed = mindist;
        closestverts = sub_borderinds(mini);
        closestverts_incluster = logical(groupmatches(closestverts,i));
        mindist_directed(closestverts_incluster) = mindist_directed(closestverts_incluster) * -1;
        border_movement_directed(group_borderinds,s) = mindist_directed;
        
        
        pct_overlap(logical(groupmatches(:,i)),s) = sum(SA_LR(logical(groupmatches(:,i) .* submatches(:,i)))) ./ sum(SA_LR(logical(groupmatches(:,i))));
        size_change(logical(groupmatches(:,i)),s) = (sum(SA_LR(logical(submatches(:,i)))) - sum(SA_LR(logical(groupmatches(:,i))))) ./ sum(SA_LR(logical(groupmatches(:,i))));
        %size_ratio(logical(groupmatches(:,i)),s) = nnz(submatches(:,i)) ./ nnz(groupmatches(:,i));
        
        [ign, groupcentroid] = min(mean(distances(:,logical(groupmatches(:,i))),2));
        [ign, subcentroid] = min(mean(distances(:,logical(submatches(:,i))),2));
        centroid_movement(logical(groupmatches(:,i)),s) = distances(groupcentroid,subcentroid);
        
        
%         hemnum = 1;
%         if all(find(groupmatches(:,i)) > nnz(mask{1}))
%             hemnum = 2;
%         end
%        
%         hemgroupcluster = zeros(size(mask{hemnum}));
%         hemgroupcluster(logical(mask{hemnum})) = groupmatches((1:nnz(mask{hemnum})) + (nnz(mask{1})*(hemnum-1)),i);
%         groupinds = find(hemgroupcluster);
%         groupcentroid = mean(verts{hemnum}(groupinds,:),1);
%         
%         hemsubcluster = zeros(size(mask{hemnum}));
%         hemsubcluster(logical(mask{hemnum})) = submatches((1:nnz(mask{hemnum})) + (nnz(mask{1})*(hemnum-1)),i);
%         subinds = find(hemsubcluster);
%         subcentroid = mean(verts{hemnum}(subinds,:),1);
%         
%         centroid_movement(logical(groupmatches(:,i)),s) = pdist([groupcentroid;subcentroid]);
        
    end
end

% out = zeros(size(mostcommon));
% out(1:nsurfverts) = nanmean(pct_overlap,2);
% cifti_write_wHDR(out,[],'Avg_Percent_Cluster_Overlap')
% 
% out = zeros(size(mostcommon));
% out(1:nsurfverts) = nanmean(centroid_movement,2);
% cifti_write_wHDR(out,[],'Avg_Cluster_Centroid_Movement')

out = zeros(size(mostcommon));
out(1:nsurfverts) = nanmean(border_movement,2);
cifti_write_wHDR(out,[],'Avg_Cluster_Border_Movement')

out = zeros(size(mostcommon));
out(1:nsurfverts) = nanmean(abs(size_change),2);
cifti_write_wHDR(out,[],'Avg_Cluster_Size_Change')

out = zeros(size(mostcommon),size(networkconnection_bysub,2));
out(1:nsurfverts,:) = size_change;
cifti_write_wHDR(out,[],'Cluster_Size_Change_bysub')

out = zeros(size(mostcommon,1),size(networkconnection_bysub,2));
out(1:nsurfverts,:) = border_movement_directed;
cifti_write_wHDR(out,[],'Border_Movement_bysub')

%save('border_distance_directed_bysub.mat','border_movement_directed','-v7.3')



%% Border Segment Movement Correlations

ncortverts = 59412;

border_movement_directed = cifti_read('Border_Movement_bysub.dtseries.nii');
border_movement_directed = border_movement_directed(1:ncortverts,:);


neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');

mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');

edges = zeros(size(mostcommon));

for i = 1:ncortverts
    
    vertneighs = neighbors(i,2:end); vertneighs(isnan(vertneighs)) = [];
    if any(mostcommon(vertneighs) == mostcommon(i))
        edges(i) = mostcommon(i);
    end
end

cifti_write_wHDR(edges,[],'Variability_LR_consensus_dice_mostcommon_clean50mm_edges')

IDs = unique(mostcommon); IDs(IDs==0) = [];

segments = zeros(size(edges,1),0);

segmentIDs = zeros(0,2);

segment_movement = zeros(0,size(border_movement_directed,2));

for ID1 = IDs(:)'
    for ID2 = IDs(:)'
        if ID1 ~= ID2
            
            ID1_nexttoID2_edges = zeros(size(edges));
            
            ID1_edgeverts = find(edges==ID1);
            
            for vert = ID1_edgeverts(:)'
                
                vertneighs = neighbors(vert,2:end); vertneighs(isnan(vertneighs)) = [];
                
                if any(edges(vertneighs)==ID2)
                    
                    ID1_nexttoID2_edges(vert) = 1;
                end
            end
            
            temp = metric_cluster_cifti(ID1_nexttoID2_edges,.5,1.5,0);
            
            segments(:,end+1:end+size(temp,2)) = temp;
            
            segmentIDs(end+1:end+size(temp,2),:) = repmat([ID1 ID2],size(temp,2),1);
            
            for segnum = 1:size(temp,2)
                segment_movement(end+1,:) = mean(border_movement_directed(logical(temp(:,segnum)),:),1);
            end
            
        end
    end
end

segment_movement_correlations = zeros(size(segment_movement,1));
segment_movement_significance = zeros(size(segment_movement,1));
for i = 1:size(segment_movement,1)
    for j = 1:size(segment_movement,1)
        inds = logical((~isnan(segment_movement(i,:))) .* (~isnan(segment_movement(j,:))));
        
        [R, P] = corrcoef(segment_movement(i,inds)',segment_movement(j,inds)');
        
        segment_movement_correlations(i,j) = R(1,2);
        segment_movement_significance(i,j) = P(1,2);
    end
end
%thresh = .05;
thresh = .05 / (length(segment_movement_correlations) / 2 *(length(segment_movement_correlations)-1));
segment_movement_correlations = segment_movement_correlations .* (segment_movement_significance < thresh);

figure
parcel_correlmat_figmaker(segment_movement_correlations,segmentIDs(:,1),[-.6 .6],'Border Movement Correlations')
export_fig(gcf,'Border Movement Correlations by Segment.tif')



%%
% 
% centroid_movement_regresssize = ones(size(centroid_movement)) * NaN;
% 
% for vert = 1:size(centroid_movement,1)
%     nonnaninds = ~isnan(centroid_movement(vert,:));
%     [B,BINT,R] = regress(centroid_movement(vert,nonnaninds)',abs(log(size_change(vert,nonnaninds)')));
%     centroid_movement_regresssize(vert,nonnaninds) = R;
% end
% 
% out = zeros(size(mostcommon));
% out(1:nsurfverts) = nanmean(centroid_movement_regresssize,2);
% cifti_write_wHDR(out,[],'Avg_Cluster_Centroid_Movement_Regressize')


%% Cluster Metric Correlations

IDs = unique(mostcommon); IDs(IDs==0) = [];
mostcommon_clusters = zeros(length(mostcommon),0);
for ID = IDs(:)'
outputcifti = metric_cluster_cifti_surfacearea(mostcommon,ID-.01,ID+.01,10);
mostcommon_clusters(:,end+1:end+size(outputcifti,2)) = outputcifti;
end

pct_overlap_byclus = zeros(size(networkconnection_bysub,2),size(mostcommon_clusters,2));
centroid_movement_byclus = zeros(size(networkconnection_bysub,2),size(mostcommon_clusters,2));
border_movement_byclus = zeros(size(networkconnection_bysub,2),size(mostcommon_clusters,2));
size_change_byclus = zeros(size(networkconnection_bysub,2),size(mostcommon_clusters,2));
for i = 1:size(mostcommon_clusters,2)
    pct_overlap_byclus(:,i) = mean(pct_overlap(logical(mostcommon_clusters(:,i)),:),1);
    centroid_movement_byclus(:,i) = mean(centroid_movement(logical(mostcommon_clusters(:,i)),:),1);
    border_movement_byclus(:,i) = nanmean(border_movement(logical(mostcommon_clusters(:,i)),:),1);
    size_change_byclus(:,i) = nanmean(size_change(logical(mostcommon_clusters(:,i)),:),1);
end

pct_overlap_rs = zeros(size(mostcommon_clusters,2),size(mostcommon_clusters,2));
pct_overlap_ps = pct_overlap_rs; centroid_movement_rs = pct_overlap_rs; centroid_movement_ps = pct_overlap_rs; border_movement_rs = pct_overlap_rs; border_movement_ps = pct_overlap_rs; size_change_rs = pct_overlap_rs; size_change_ps = pct_overlap_rs;

count_matched = zeros(size(mostcommon));

for i = 1:size(mostcommon_clusters,2)
    nonnaninds_i = ~isnan(pct_overlap_byclus(:,i));
    count_matched(logical(mostcommon_clusters(:,i))) = nnz(nonnaninds_i);
    for j = 1:size(mostcommon_clusters,2)
        if i ~= j
        nonnaninds_j = ~isnan(pct_overlap_byclus(:,j));
        nonnaninds_ij = logical(nonnaninds_i .* nonnaninds_j);
        if nnz(nonnaninds_ij) > 2
            [pct_overlap_rs(i,j), pct_overlap_ps(i,j)] = corr(pct_overlap_byclus(nonnaninds_ij,i),pct_overlap_byclus(nonnaninds_ij,j));
            [centroid_movement_rs(i,j), centroid_movement_ps(i,j)] = corr(centroid_movement_byclus(nonnaninds_ij,i),centroid_movement_byclus(nonnaninds_ij,j));
            [border_movement_rs(i,j), border_movement_ps(i,j)] = corr(border_movement_byclus(nonnaninds_ij,i),border_movement_byclus(nonnaninds_ij,j));
            [size_change_rs(i,j), size_change_ps(i,j)] = corr(size_change_byclus(nonnaninds_ij,i),size_change_byclus(nonnaninds_ij,j));
        end
        end
    end
end

thresh = .05 / nnz(tril(ones(size(pct_overlap_rs)),1));

%cifti_write_wHDR(count_matched,[],'Matched_Cluster_Count')

% disp('Pct Overlap')
% disp(pct_overlap_rs .* (pct_overlap_ps < thresh))
% 
% disp('Centroid Movement')
% disp(centroid_movement_rs .* (centroid_movement_ps < thresh))
% 
% disp('Border Movement')
% disp(border_movement_rs .* (border_movement_ps < thresh))
% 
% disp('Size Ratio')
% disp(size_ratio_rs .* (size_ratio_ps < thresh))

assignments = zeros(size(mostcommon_clusters,2),1);
for i = 1:length(assignments)
    assignments(i) = mean(mostcommon(logical(mostcommon_clusters(:,i))));
end

figure
parcel_correlmat_figmaker((pct_overlap_rs .* (pct_overlap_ps < thresh)),assignments,[-.6 .6],'Percent Overlap Correlations')
export_fig(gcf,'Percent Overlap Correlations.tif')

figure
parcel_correlmat_figmaker((centroid_movement_rs .* (centroid_movement_ps < thresh)),assignments,[-.6 .6],'Centroid Movement Correlations')
export_fig(gcf,'Centroid Movement Correlations.tif')

figure
parcel_correlmat_figmaker((border_movement_rs .* (border_movement_ps < thresh)),assignments,[-.6 .6],'Border Movement Correlations')
export_fig(gcf,'Border Movement Correlations.tif')

figure
parcel_correlmat_figmaker((size_change_rs .* (size_change_ps < thresh)),assignments,[-.6 .6],'Size Change Correlations')
export_fig(gcf,'Size Change Correlations.tif')

  
%% Cluster Size and Border Movement Regressing data quality

ncortverts = 59412;

mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');
IDs = unique(mostcommon); IDs(IDs==0) = [];
mostcommon_clusters = zeros(length(mostcommon),0);
for ID = IDs(:)'
outputcifti = metric_cluster_cifti_surfacearea(mostcommon,ID-.01,ID+.01,10);
mostcommon_clusters(:,end+1:end+size(outputcifti,2)) = outputcifti;
end

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

sulcdepth_diff = zeros(ncortverts,length(subjects));
sulcdepth_diff_clus = zeros(size(mostcommon_clusters,2),length(subjects));
curv_diff = zeros(ncortverts,length(subjects));
curv_diff_clus = zeros(size(mostcommon_clusters,2),length(subjects));
ntimepoints = zeros(length(subjects),1);
meanFD = zeros(length(subjects),1);

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
sulcdepth_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_L = sulcdepth_group_L.cdata(logical(maskL));
curv_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.curvature.32k_fs_LR.shape.gii'); curv_group_L = curv_group_L.cdata(logical(maskL));

maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
sulcdepth_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_R = sulcdepth_group_R.cdata(logical(maskR));
curv_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.curvature.32k_fs_LR.shape.gii'); curv_group_R = curv_group_R.cdata(logical(maskR));


for s = 1:length(subjects)
    sulcdepth_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_diff(:,s) = [sulcdepth_L.cdata(logical(maskL)) - sulcdepth_group_L ; sulcdepth_R.cdata(logical(maskR)) - sulcdepth_group_R];
        
    curv_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.curvature.32k_fs_LR.shape.gii']);
    curv_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.curvature.32k_fs_LR.shape.gii']);
    curv_diff(:,s) = [curv_L.cdata(logical(maskL)) - curv_group_L ; curv_R.cdata(logical(maskR)) - curv_group_R];
    
    for clusnum = 1:size(mostcommon_clusters,2)
        sulcdepth_diff_clus(clusnum,s) = mean(sulcdepth_diff(logical(mostcommon_clusters(:,clusnum)),s));
        curv_diff_clus(clusnum,s) = mean(curv_diff(logical(mostcommon_clusters(:,clusnum)),s));
    end
    
    tmask = load(tmasks{s});
    ntimepoints(s) = nnz(tmask);
    
    fd_file = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/FCPROCESS_bandpass_interp_nosmooth/' subjects{s} '/total_FD.txt'];
    if exist(fd_file)
        FD = load(fd_file);
        meanFD(s) = mean(FD(logical(tmask)));
    else
        fd_file = ['/data/cn5/selfRegulation/V4Process/initial_process/' subjects{s} '/' subjects{s} '/total_FD.txt'];
        FD = load(fd_file);
        meanFD(s) = mean(FD(logical(tmask)));
    end
    
end

size_change = abs(cifti_read('Cluster_Size_Change_bysub.dtseries.nii'));

Size_V_Sulc_Curv_Ntimepoints_FD_beta = zeros(66697,4);
Size_V_Sulc_Curv_Ntimepoints_FD_pval = zeros(66697,4);
Size_V_Sulc_Curv_Ntimepoints_FD_resid = ones(66697,size(size_change,2)) * NaN;

for vert = 1:ncortverts
    
    clusnum_ind = find(mostcommon_clusters(vert,:));
    
    inds = logical(~isnan(size_change(vert,:)));
    
    STATS = regstats(size_change(vert,inds)',[sulcdepth_diff_clus(clusnum_ind,inds)' curv_diff_clus(clusnum_ind,inds)' ntimepoints(inds) meanFD(inds)],'linear',{'tstat','r'});
    Size_V_Sulc_Curv_Ntimepoints_FD_beta(vert,:) = STATS.tstat.beta(2:end);
    Size_V_Sulc_Curv_Ntimepoints_FD_pval(vert,:) = STATS.tstat.pval(2:end);
    Size_V_Sulc_Curv_Ntimepoints_FD_resid(vert,inds) = STATS.r;
    
end

cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_beta,[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_beta')
cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_pval,[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_pval')
cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_resid,[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_resid')
cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_resid .* sign(cifti_read('Cluster_Size_Change_bysub.dtseries.nii')),[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_resid_directed')



border_shift = abs(cifti_read('Border_Movement_bysub.dtseries.nii'));

Border_V_Sulc_Curv_Ntimepoints_FD_beta = zeros(66697,4);
Border_V_Sulc_Curv_Ntimepoints_FD_pval = zeros(66697,4);
Border_V_Sulc_Curv_Ntimepoints_FD_resid = ones(66697,size(border_shift,2)) * NaN;

verts = find(any(~isnan(border_shift),2) .* any(border_shift > 0 , 2));

for vert = verts(:)'
    
    inds = logical(~isnan(border_shift(vert,:)));
    
    STATS = regstats(border_shift(vert,inds)',[sulcdepth_diff(vert,inds)' curv_diff(vert,inds)' ntimepoints(inds) meanFD(inds)],'linear',{'tstat','r'});
    Border_V_Sulc_Curv_Ntimepoints_FD_beta(vert,:) = STATS.tstat.beta(2:end);
    Border_V_Sulc_Curv_Ntimepoints_FD_pval(vert,:) = STATS.tstat.pval(2:end);
    Border_V_Sulc_Curv_Ntimepoints_FD_resid(vert,inds) = STATS.r;
    
end

cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_beta,[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_beta')
cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_pval,[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_pval')
cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_resid,[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid')
cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_resid .* sign(cifti_read('Border_Movement_bysub.dtseries.nii')),[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid_directed')



%% Residual Cluster Size and Border Movement correlations

ncortverts = 59412;

mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');
IDs = unique(mostcommon); IDs(IDs==0) = [];
mostcommon_clusters = zeros(length(mostcommon),0);
for ID = IDs(:)'
outputcifti = metric_cluster_cifti_surfacearea(mostcommon,ID-.01,ID+.01,10);
mostcommon_clusters(:,end+1:end+size(outputcifti,2)) = outputcifti;
end

border_movement_resid = cifti_read('BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid_directed.dtseries.nii');
size_change_resid = cifti_read('Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_resid_directed.dtseries.nii');

border_movement_resid_byclus = zeros(size(border_movement_resid,2),size(mostcommon_clusters,2));
size_change_resid_byclus = zeros(size(size_change_resid,2),size(mostcommon_clusters,2));

for i = 1:size(mostcommon_clusters,2)
    border_movement_resid_byclus(:,i) = nanmean(border_movement_resid(logical(mostcommon_clusters(:,i)),:),1);
    size_change_resid_byclus(:,i) = nanmean(size_change_resid(logical(mostcommon_clusters(:,i)),:),1);
end

border_movement_resid_rs = zeros(size(mostcommon_clusters,2),size(mostcommon_clusters,2));
border_movement_resid_ps = border_movement_resid_rs; size_change_resid_rs = border_movement_resid_rs; size_change_resid_ps = border_movement_resid_rs;

count_matched = zeros(size(mostcommon));

for i = 1:size(mostcommon_clusters,2)
    nonnaninds_i = ~isnan(border_movement_resid_byclus(:,i));
    count_matched(logical(mostcommon_clusters(:,i))) = nnz(nonnaninds_i);
    for j = 1:size(mostcommon_clusters,2)
        if i ~= j
        nonnaninds_j = ~isnan(border_movement_resid_byclus(:,j));
        nonnaninds_ij = logical(nonnaninds_i .* nonnaninds_j);
        if nnz(nonnaninds_ij) > 2
            [border_movement_resid_rs(i,j), border_movement_resid_ps(i,j)] = corr(border_movement_resid_byclus(nonnaninds_ij,i),border_movement_resid_byclus(nonnaninds_ij,j));
            [size_change_resid_rs(i,j), size_change_resid_ps(i,j)] = corr(size_change_resid_byclus(nonnaninds_ij,i),size_change_resid_byclus(nonnaninds_ij,j));
        end
        end
    end
end

thresh = .05 / nnz(tril(ones(size(border_movement_resid_rs)),1));


assignments = zeros(size(mostcommon_clusters,2),1);
for i = 1:length(assignments)
    assignments(i) = mean(mostcommon(logical(mostcommon_clusters(:,i))));
end

figure
parcel_correlmat_figmaker((border_movement_resid_rs .* (border_movement_resid_ps < thresh)),assignments,[-.6 .6],'Residual Border Movement Correlations')
export_fig(gcf,'Residual Border Movement Correlations.tif')

figure
parcel_correlmat_figmaker((size_change_resid_rs .* (size_change_resid_ps < thresh)),assignments,[-.6 .6],'Residual Size Change Correlations')
export_fig(gcf,'Residual Size Change Correlations.tif')

%% Residual Border Segment movement correlations

ncortverts = 59412;

border_movement_directed = cifti_read('BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid_directed.dtseries.nii');
border_movement_directed = border_movement_directed(1:ncortverts,:);


neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');

mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');

edges = zeros(size(mostcommon));

for i = 1:ncortverts
    
    vertneighs = neighbors(i,2:end); vertneighs(isnan(vertneighs)) = [];
    if any(mostcommon(vertneighs) == mostcommon(i))
        edges(i) = mostcommon(i);
    end
end

cifti_write_wHDR(edges,[],'Variability_LR_consensus_dice_mostcommon_clean50mm_edges')

IDs = unique(mostcommon); IDs(IDs==0) = [];

segments = zeros(size(edges,1),0);

segmentIDs = zeros(0,2);

segment_movement = zeros(0,size(border_movement_directed,2));

for ID1 = IDs(:)'
    for ID2 = IDs(:)'
        if ID1 ~= ID2
            
            ID1_nexttoID2_edges = zeros(size(edges));
            
            ID1_edgeverts = find(edges==ID1);
            
            for vert = ID1_edgeverts(:)'
                
                vertneighs = neighbors(vert,2:end); vertneighs(isnan(vertneighs)) = [];
                
                if any(edges(vertneighs)==ID2)
                    
                    ID1_nexttoID2_edges(vert) = 1;
                end
            end
            
            temp = metric_cluster_cifti(ID1_nexttoID2_edges,.5,1.5,0);
            
            segments(:,end+1:end+size(temp,2)) = temp;
            
            segmentIDs(end+1:end+size(temp,2),:) = repmat([ID1 ID2],size(temp,2),1);
            
            for segnum = 1:size(temp,2)
                segment_movement(end+1,:) = mean(border_movement_directed(logical(temp(:,segnum)),:),1);
            end
            
        end
    end
end

segment_movement_correlations = zeros(size(segment_movement,1));
segment_movement_significance = zeros(size(segment_movement,1));
for i = 1:size(segment_movement,1)
    for j = 1:size(segment_movement,1)
        inds = logical((~isnan(segment_movement(i,:))) .* (~isnan(segment_movement(j,:))));
        
        [R, P] = corrcoef(segment_movement(i,inds)',segment_movement(j,inds)');
        
        segment_movement_correlations(i,j) = R(1,2);
        segment_movement_significance(i,j) = P(1,2);
    end
end
%thresh = .05;
thresh = .05 / (length(segment_movement_correlations) / 2 *(length(segment_movement_correlations)-1));
segment_movement_correlations = segment_movement_correlations .* (segment_movement_significance < thresh);

figure
parcel_correlmat_figmaker(segment_movement_correlations,segmentIDs(:,1),[-.6 .6],'Residual Border Movement Correlations')
export_fig(gcf,'Residual Border Movement Correlations by Segment.tif')


%% Variable Region Correlations

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
alternate_ID_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');
mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');

alternate_ID_pct_bysub_byregion = zeros(size(networkconnection_bysub,2),size(alternate_ID_regions,2));
assignments = zeros(size(alternate_ID_regions,2),1);

for r = 1:size(alternate_ID_pct_bysub_byregion,2)
    assignments(r) = mode(mostcommon(logical(alternate_ID_regions(:,r))));
    for s = 1:size(alternate_ID_pct_bysub_byregion,1)
        alternate_ID_pct_bysub_byregion(s,r) = nnz(networkconnection_bysub(logical(logical(alternate_ID_regions(:,r))),s)==max(alternate_ID_regions(:,r))) / nnz(alternate_ID_regions(:,r));
    end
end

thresh = .05 / (size(alternate_ID_pct_bysub_byregion,2) / 2 *(size(alternate_ID_pct_bysub_byregion,2)-1));

[alternate_ID_region_correlations, alternate_ID_region_significance] = corrcoef(alternate_ID_pct_bysub_byregion);

alternate_ID_region_correlations_thresh = alternate_ID_region_correlations .* (alternate_ID_region_significance < thresh);

figure
parcel_correlmat_figmaker(alternate_ID_region_correlations_thresh,assignments,[-.6 .6],'Variable Region Correlations')


%% Variable Region Correlations Confound Regressed

ncortverts = 59412;

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
alternate_ID_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');
mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');

alternate_ID_pct_bysub_byregion = zeros(size(networkconnection_bysub,2),size(alternate_ID_regions,2));
alternate_ID_pct_bysub_byregion_resid = zeros(size(networkconnection_bysub,2),size(alternate_ID_regions,2));

alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_beta = zeros(size(mostcommon,1),4);
alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_pval = zeros(size(mostcommon,1),4);

assignments = zeros(size(alternate_ID_regions,2),1);

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

sulcdepth_diff = zeros(ncortverts,length(subjects));
curv_diff = zeros(ncortverts,length(subjects));
ntimepoints = zeros(length(subjects),1);
meanFD = zeros(length(subjects),1);

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
sulcdepth_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_L = sulcdepth_group_L.cdata(logical(maskL));
curv_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.curvature.32k_fs_LR.shape.gii'); curv_group_L = curv_group_L.cdata(logical(maskL));

maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
sulcdepth_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_R = sulcdepth_group_R.cdata(logical(maskR));
curv_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.curvature.32k_fs_LR.shape.gii'); curv_group_R = curv_group_R.cdata(logical(maskR));

for s = 1:size(alternate_ID_pct_bysub_byregion,1)
    sulcdepth_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_diff(:,s) = [sulcdepth_L.cdata(logical(maskL)) - sulcdepth_group_L ; sulcdepth_R.cdata(logical(maskR)) - sulcdepth_group_R];
        
    curv_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.curvature.32k_fs_LR.shape.gii']);
    curv_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.curvature.32k_fs_LR.shape.gii']);
    curv_diff(:,s) = [curv_L.cdata(logical(maskL)) - curv_group_L ; curv_R.cdata(logical(maskR)) - curv_group_R];
    
    
    
    tmask = load(tmasks{s});
    ntimepoints(s) = nnz(tmask);
    
    fd_file = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/FCPROCESS_bandpass_interp_nosmooth/' subjects{s} '/total_FD.txt'];
    if exist(fd_file)
        FD = load(fd_file);
        meanFD(s) = mean(FD(logical(tmask)));
    else
        fd_file = ['/data/cn5/selfRegulation/V4Process/initial_process/' subjects{s} '/' subjects{s} '/total_FD.txt'];
        FD = load(fd_file);
        meanFD(s) = mean(FD(logical(tmask)));
    end
end


for r = 1:size(alternate_ID_pct_bysub_byregion,2)
    
    regioninds = logical(alternate_ID_regions(:,r));
    
    assignments(r) = mode(mostcommon(regioninds));
    for s = 1:size(alternate_ID_pct_bysub_byregion,1)
        alternate_ID_pct_bysub_byregion(s,r) = nnz(networkconnection_bysub(logical(regioninds),s)==max(alternate_ID_regions(:,r))) / nnz(regioninds);
    end
    
    STATS = regstats(alternate_ID_pct_bysub_byregion(:,r)',[mean(sulcdepth_diff(regioninds,:),1)' mean(curv_diff(regioninds,:),1)' ntimepoints meanFD],'linear',{'tstat','r'});
    alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_beta(regioninds,:) = repmat(STATS.tstat.beta(2:end)',nnz(regioninds),1);
    alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_pval(regioninds,:) = repmat(STATS.tstat.pval(2:end)',nnz(regioninds),1);
    alternate_ID_pct_bysub_byregion_resid(:,r) = STATS.r;
    
end

thresh = .05 / (size(alternate_ID_pct_bysub_byregion,2) / 2 *(size(alternate_ID_pct_bysub_byregion,2)-1));

[alternate_ID_region_correlations, alternate_ID_region_significance] = corrcoef(alternate_ID_pct_bysub_byregion_resid);

alternate_ID_region_correlations_thresh = alternate_ID_region_correlations .* (alternate_ID_region_significance < thresh);

figure
parcel_correlmat_figmaker(alternate_ID_region_correlations_thresh,assignments,[-.6 .6],'Residual Variable Region Correlations')

cifti_write_wHDR(alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_beta,[],'Variable_Regions_V_Sulc_Curv_Ntimepoints_FD_beta')
cifti_write_wHDR(alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_pval,[],'Variable_Regions_V_Sulc_Curv_Ntimepoints_FD_pval')
    

%% Variable region connectivity t-tests


percentregion_thresh = .5;

all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8_selected.dtseries.nii');
variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon.dtseries.nii');

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

ncortverts = 59412;

%main = zeros(ncortverts,size(variable_regions,2));
%alternate = zeros(ncortverts,size(variable_regions,2));
alternatecount = zeros(1,size(variable_regions,2));

for s = 1:length(subjects)
    disp(subjects{s})
    tmask = load(tmasks{s});
    subdata = cifti_read(ciftifiles{s});
    subdata = subdata(:,logical(tmask));
    
    for r = 1:size(variable_regions,2)
        
        if s==1
            alternate{r} = zeros(ncortverts,0);
            main{r} = zeros(ncortverts,0);
        end
        
        mainID = mode(mostcommon(logical(variable_regions(:,r)))); 
        altID = mean(all_variable_regions(logical(variable_regions(:,r))));
        if (nnz(networkconnection_bysub(logical(variable_regions(:,r)),s)==altID) / nnz(variable_regions(:,r))) >= percentregion_thresh
            alternatecount(r) = alternatecount(r)+1;
            
            indices = logical(variable_regions(:,r));% .* (networkconnection_bysub(:,s)==altID));
            
            corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
            corrpattern(isnan(corrpattern)) = 0;
            
            %alternate(:,r) = alternate(:,r) + corrpattern(1:ncortverts);
            alternate{r}(:,end+1) = corrpattern(1:ncortverts);
            
        else
            
            indices = logical(variable_regions(:,r));% .* (networkconnection_bysub(:,s)==mainID));
            
            corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
            corrpattern(isnan(corrpattern)) = 0;
            
            %main(:,r) = main(:,r) + corrpattern(1:ncortverts);
            main{r}(:,end+1) = corrpattern(1:ncortverts);
        end
    end
end

mainout = zeros(size(mostcommon,1),size(variable_regions,2));
altout = zeros(size(mostcommon,1),size(variable_regions,2));
tout = zeros(size(mostcommon,1),size(variable_regions,2));
tthresh = tinv(1-(.05 / ncortverts / size(variable_regions,2)),length(subjects));
for r = 1:size(variable_regions,2)
    
    %     alternate(:,r) = alternate(:,r) ./ alternatecount(r);
    %     main(:,r) = main(:,r) ./ (length(subjects) - alternatecount(r));
    
    mainout(1:ncortverts,r) = mean(main{r},2);
    altout(1:ncortverts,r) = mean(alternate{r},2);
    
    [H,P,CI,STATS] = ttest2(main{r}',alternate{r}');
    tout(1:ncortverts,r) = STATS.tstat;
    
end

disp(alternatecount)

% cifti_write_wHDR(mainout,[],'MeanConnectivity_mainID_byregion')
% cifti_write_wHDR(altout,[],'MeanConnectivity_alternateID_byregion')
% cifti_write_wHDR(tout,[],['MeanConnectivity_T_MainVsAlternateID_byregion_correctedT' num2str(tthresh)])
% cifti_write_wHDR((mainout - altout),[],'MeanConnectivity_Diff_MainVsAlternateID_byregion')

cifti_write_wHDR(mainout,[],'MeanConnectivity_mainID_bywholeregion')
cifti_write_wHDR(altout,[],'MeanConnectivity_alternateID_bywholeregion')
cifti_write_wHDR(tout,[],['MeanConnectivity_T_MainVsAlternateID_bywholeregion_correctedT' num2str(tthresh)])
cifti_write_wHDR((mainout - altout),[],'MeanConnectivity_Diff_MainVsAlternateID_bywholeregion')


%% Variable region characteristics

distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
SA = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas.dtseries.nii');
coords = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_ASCIIformat_coords.dtseries.nii');

Variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');

for i = 1:size(Variable_regions,2)
    inds = find(Variable_regions(:,i));
    disp(['Region ' num2str(i)])
    this_SA = sum(SA(inds));
    disp(['Surface area: ' num2str(this_SA)])

    inddistances = sum(distances(inds,inds),1);
    [ign mini] = min(inddistances);
    centroid = inds(mini);
    disp(['Centroid: ' num2str(coords(centroid,:))])
end

