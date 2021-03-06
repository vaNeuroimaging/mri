%% Probabilistic system maps
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
IDs = unique(networkconnection_bysub); IDs(IDs==0) = [];

prob_maps = zeros(size(networkconnection_bysub,1),length(IDs));

for IDnum = 1:length(IDs)
    ID = IDs(IDnum);
    this_prob_map = sum((networkconnection_bysub==ID),2) ./ size(networkconnection_bysub,2);
    prob_maps(:,IDnum) = this_prob_map;
end

cifti_write_wHDR(prob_maps,[],'Probabilistic_system_maps')



%% Alternate Connected Regions Close and Distant to System Borders
subthresh = 46;
distance = 8;
sizethreshmm = 100;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

hems = {'L','R'};
for hemnum = 1:length(hems)
    hem = hems{hemnum};
    %edges = gifti(['Variability_' hem '_consensus_dice_mostcommon_edges.func.gii']); edges = edges.cdata;
    edges = gifti(['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_' hem '_edges.func.gii']); edges = edges.cdata;
    variability = gifti(['Variability_' hem '_consensus_dice.func.gii']); variability = variability.cdata(:,2:end);
    numsubs = gifti(['Variability_' hem '_consensus_dice_clustersize.func.gii']); numsubs = numsubs.cdata(:,2:end);
    variability(numsubs<subthresh) = 0;
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])
    
    distvariableregions = zeros(32492,1);
    closevariableregions = ones(32492,1);
    variabilityfromedge = ones(32491,1) .* -1;
    
    IDs = unique(variability); IDs(IDs==0) = [];
    for ID = IDs(:)'
        disp(ID)
        thisID_distvariableregions = zeros(32492,1);
        
        inds = [];
        for col = 1:size(variability,2)
            inds = [inds ;find(variability(:,col)==ID)];
        end
        
        for ind = inds(:)'
            withindist = geo_distances(:,ind) <= distance;
            if ~any(edges(withindist)==ID)
                thisID_distvariableregions(ind) = 1;
            end
        end
        
        
        
        temp = metric_cluster_surfacearea(thisID_distvariableregions,0.5,1.5,sizethreshmm,hem);
        for i = 1:size(temp,2)
            distvariableregions(logical(temp(:,i))) = ID;
        end
        
        
        closevariableregions(logical(thisID_distvariableregions)) = 0;
        
        thisIDedgeinds = find(edges==ID);
        thisIDvariableverts = any(variability==ID,2);
        thisIDvariableverts(logical(thisID_distvariableregions)) = 0;
        thisIDvariableverts_clusters = metric_cluster_surfacearea(thisIDvariableverts,.5,1.5,10,hem);
        clusterborders = zeros(size(thisIDvariableverts_clusters));
        for clusternum = 1:size(thisIDvariableverts_clusters,2)
            
            clusterinds = find(thisIDvariableverts_clusters(:,clusternum)==0);
            for ind = clusterinds(:)'
                indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
                if any(thisIDvariableverts_clusters(indneighs,clusternum)==1) && ~(edges(ind)==ID) && ~(any(edges(indneighs)==ID))
                    clusterborders(ind,clusternum) = 1;
                end
            end
        end
        
        for edgeind = thisIDedgeinds(:)'
            adjacentcluster = 0;
            indneighs = neighbors(edgeind,2:7); indneighs(isnan(indneighs)) = [];
            for clusternum = 1:size(thisIDvariableverts_clusters,2)
                if any(intersect(indneighs,find(thisIDvariableverts_clusters(:,clusternum))))
                    adjacentcluster = clusternum;
                    break
                end
            end
            if adjacentcluster > 0  && (nnz(clusterborders(:,adjacentcluster)) > 0)
                mindist = min(geo_distances(edgeind,logical(clusterborders(:,adjacentcluster))));
                variabilityfromedge(edgeind) = mindist;
            end
        end
            
            
                
        
    end
    
    save(gifti(single(distvariableregions)),['Variabile_regions_' hem '_consensus_dice_distance' num2str(distance) '.func.gii'])
    %save(gifti(single(variabilityfromedge)),['Variabile_regions_' hem '_distance_from_edge.func.gii'])
    
        
    
end

gifti_to_cifti(['Variabile_regions_L_consensus_dice_distance' num2str(distance) '.func.gii'],['Variabile_regions_R_consensus_dice_distance' num2str(distance) '.func.gii'],['Variabile_regions_LR_consensus_dice_distance' num2str(distance)])
colored = cifti_read(['Variabile_regions_LR_consensus_dice_distance' num2str(distance) '.dtseries.nii']);
IDs = unique(colored); IDs(IDs==0) = [];
out = zeros(size(colored,1),0);
for ID = IDs(:)'
    temp = metric_cluster_cifti(colored,ID-.5,ID+.5,0);
    out(:,end+1:end+size(temp,2)) = temp;
end

out = out .* repmat(colored,1,size(out,2));
cifti_write_wHDR(out,[],['Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated']);



%% Make Stripey Figs of Selected Regions

subthresh = 46;
distance = 8;
% selected = [2 5 9 10 13 14 18 22 26 29 30 31 35 39 50];
% allregions = cifti_read(['Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated.dtseries.nii']);
% selectedregions = allregions(:,selected);
% cifti_write_wHDR(selectedregions,[],['Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated_selected']);
cifti_to_gifti(['Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated_selected.dtseries.nii']);


hems = {'L','R'};
for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    selectedregions_hem = gifti(['Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated_selected_' hem '.func.gii']);
    mask = logical(sum(selectedregions_hem.cdata,2));


    save(gifti(single(mask)),['Variable_and_distant_mask_' hem '.func.gii'])
    system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variable_and_distant_mask_' hem '.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Variable_and_distant_mask_' hem '_164.func.gii -largest'])
    distmask164 = gifti(['Variable_and_distant_mask_' hem '_164.func.gii']); distmask164 = distmask164.cdata;
    
    
    %save(gifti(single(logical(closevariableregions))),['Variable_and_close_mask_' hem '.func.gii'])
    %system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variable_and_close_mask_' hem '.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Variable_and_close_mask_' hem '_164.func.gii -largest'])
    %closemask164 = gifti(['Variable_and_close_mask_' hem '_164.func.gii']); closemask164 = closemask164.cdata;
    
    %variability = gifti(['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned_' hem '.func.gii']);
    %save(gifti(single(variability.cdata(:,1))),['Variability_' hem '_consensus_dice_mostcommon.func.gii'])
    %system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variability_' hem '_consensus_dice_mostcommon.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Variability_' hem '_consensus_dice_mostcommon_164.func.gii -largest'])
    nonvariable = gifti(['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned_' hem '_164.func.gii']); nonvariable = nonvariable.cdata;
    %nonvariable = gifti(['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_' hem '_164.func.gii']); nonvariable = nonvariable.cdata;
    
    variability164 = gifti(['Variability_' hem '_consensus_diceminsize_' num2str(subthresh) '_combineclusters_164.func.gii']); variability164 = variability164.cdata;
    
    distantvariable = nonvariable; distantvariable(logical(distmask164)) = variability164(logical(distmask164));
    save(gifti(single(distantvariable)),['Variability_' hem '_distantonly_minsize_' num2str(subthresh) '_164.func.gii'])
    
    system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated_selected_' hem '.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated_selected_' hem '_164.func.gii -largest'])
    
    separated_distant = gifti(['Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated_selected_' hem '_164.func.gii']); separated_distant = separated_distant.cdata;
    separated_distant_variable = logical(separated_distant) .* repmat(variability164,1,size(separated_distant,2));
    save(gifti(single(separated_distant_variable)),['Variability_' hem '_distantonly_minsize_' num2str(subthresh) '_separated_164.func.gii'])
    
%     distantonlyvariable = zeros(size(distantvariable));
%     distantonlyvariable(logical(distmask164)) = variability164(logical(distmask164));
%     save(gifti(single(distantonlyvariable)),['Variability_' hem '_distantonly_nomean_minsize_' num2str(subthresh) '_164.func.gii'])
    
    %closevariable = nonvariable; closevariable(logical(closemask164)) = variability164(logical(closemask164));
    %save(gifti(single(closevariable)),['Variability_' hem '_closeonly_minsize_' num2str(subthresh) '_164.func.gii'])
end


%% Connectivity matrices of matched patches

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
groupavgsystems = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');%'Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');

orig_avgsystems = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR.dtseries.nii');
networkconnection_bysub(logical((orig_avgsystems==13) + (orig_avgsystems==14)),:) = 0;

mask{1} = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'); mask{1} = ~mask{1}.cdata;
mask{2} = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'); mask{2} = ~mask{2}.cdata;

% verts{1} = gifti('/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii'); verts{1} = verts{1}.vertices;
% verts{2} = gifti('/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii'); verts{2} = verts{2}.vertices;


nsurfverts = nnz(mask{1}) + nnz(mask{2});

neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');
distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat');
distances = distances(1:nsurfverts,1:nsurfverts);

cluster_SA_thresh = 100;

IDs = unique(groupavgsystems); IDs(IDs==0) = [];
groupavgsystems_clusters = zeros(length(groupavgsystems),0);
for ID = IDs(:)'
    outputcifti = metric_cluster_cifti_surfacearea(groupavgsystems,ID-.01,ID+.01,cluster_SA_thresh);
    groupavgsystems_clusters(:,end+1:end+size(outputcifti,2)) = outputcifti;
end

groupavgsystems_clusters = groupavgsystems_clusters(1:nsurfverts,:);

sub_patch_correls = ones(size(groupavgsystems_clusters,2),size(groupavgsystems_clusters,2),length(subjects)) .* NaN;




for s = 1:size(networkconnection_bysub,2)
    disp(['Subject ' num2str(s)])
    fprintf('   ')
    [groupmatches, submatches] = networkClusters_matchtogroup(networkconnection_bysub(:,s),distances,groupavgsystems);
    
    groupmatched_vec = zeros(1,size(groupavgsystems_clusters,2));
    for patchnum = 1:size(groupmatches,2)
        matchind = find(any(groupavgsystems_clusters .* repmat(groupmatches(:,patchnum),1,size(groupavgsystems_clusters,2)),1));
        groupmatched_vec(matchind) = patchnum;
    end
    submatched_vec = groupmatched_vec(logical(groupmatched_vec));
    groupmatched_inds = find(groupmatched_vec);
    
    groupavg_patch_correls_thissub = zeros(size(groupavgsystems_clusters,2),size(groupavgsystems_clusters,2));
    
    subdata = ft_read_cifti_mod(ciftifiles{s});
    tmask = load(tmasks{s});
    subdata = subdata.data(:,logical(tmask));
    
    patch_timecourses = zeros(length(groupmatched_inds),size(subdata,2));
    for i = 1:length(groupmatched_inds)
        patch_timecourses(i,:) = mean(subdata(logical(submatches(:,submatched_vec(i))),:),1);
    end
    patch_corrmat = paircorr_mod(patch_timecourses');
    patch_corrmat(isnan(patch_corrmat)) = 0;
    patch_corrmat = FisherTransform(patch_corrmat);
    
    sub_patch_correls(logical(groupmatched_vec),logical(groupmatched_vec),s) = patch_corrmat;
end

assignments = zeros(size(groupavgsystems_clusters,2),1);
for i = 1:length(assignments)
    assignments(i) = mean(groupavgsystems(logical(groupavgsystems_clusters(:,i))));
end

avg_sub_patch_correls = nanmean(sub_patch_correls,3);
parcel_correlmat_figmaker(avg_sub_patch_correls,assignments,[-.6 .6],'Matched patch connectiity')
    



%% Variable Region Correlations

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
alternate_ID_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');

alternate_ID_pct_bysub_byregion = zeros(size(networkconnection_bysub,2),size(alternate_ID_regions,2));
assignments = zeros(size(alternate_ID_regions,2),1);

for r = 1:size(alternate_ID_pct_bysub_byregion,2)
    %assignments(r) = mode(mostcommon(logical(alternate_ID_regions(:,r))));
    %assignments(r) = mode(alternate_ID_regions(logical(alternate_ID_regions(:,r)),r));
    for s = 1:size(alternate_ID_pct_bysub_byregion,1)
        alternate_ID_pct_bysub_byregion(s,r) = nnz(networkconnection_bysub(logical(logical(alternate_ID_regions(:,r))),s)==max(alternate_ID_regions(:,r))) / nnz(alternate_ID_regions(:,r));
    end
end

assignments = [7 7 5 1 1 1 2 2 3 3 3 1 1 3];

thresh = .05 / (size(alternate_ID_pct_bysub_byregion,2) / 2 *(size(alternate_ID_pct_bysub_byregion,2)-1));

[alternate_ID_region_correlations, alternate_ID_region_significance] = corrcoef(alternate_ID_pct_bysub_byregion);

alternate_ID_region_correlations_thresh = alternate_ID_region_correlations .* (alternate_ID_region_significance < thresh);

figure
parcel_correlmat_figmaker(alternate_ID_region_correlations_thresh,assignments,[-.6 .6],'Variable Region Correlations')


%% Variable Region Correlations Confound Regressed

ncortverts = 59412;

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
alternate_ID_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');

reorder = [4 5 6 12 13 7 8 9 10 11 14 3 1 2];
alternate_ID_regions = alternate_ID_regions(:,reorder);

mostcommon = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');

alternate_ID_pct_bysub_byregion = zeros(size(networkconnection_bysub,2),size(alternate_ID_regions,2));
alternate_ID_pct_bysub_byregion_resid = zeros(size(networkconnection_bysub,2),size(alternate_ID_regions,2));


alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta = zeros(size(mostcommon,1),6);
alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval = zeros(size(mostcommon,1),6);
alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq = zeros(size(mostcommon,1),1);

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
    
    
    AD_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.ArealDistortion_32k_fs_LR.shape.gii']);
    AD_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.ArealDistortion_32k_fs_LR.shape.gii']);
    AD(:,s) = abs([AD_L.cdata(logical(maskL)) ; AD_R.cdata(logical(maskR))]);
    
    subSNR = ft_read_cifti_mod(['/data/cn4/evan/RestingState/Ind_variability/120_108_SNR/' subjects{s} '_SNR.dtseries.nii']);
    SNR(:,s) = subSNR.data(1:59412,:);
    
    
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
    
    STATS = regstats(alternate_ID_pct_bysub_byregion(:,r)',[mean(sulcdepth_diff(regioninds,:),1)' mean(curv_diff(regioninds,:),1)' ntimepoints meanFD mean(AD(regioninds,:),1)' mean(SNR(regioninds,:),1)'],'linear',{'tstat','r','adjrsquare'});
    alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta(regioninds,:) = repmat(STATS.tstat.beta(2:end)',nnz(regioninds),1);
    alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval(regioninds,:) = repmat(STATS.tstat.pval(2:end)',nnz(regioninds),1);
    alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq(regioninds,:) = repmat(STATS.adjrsquare,nnz(regioninds),1);
    alternate_ID_pct_bysub_byregion_resid(:,r) = STATS.r;
    
    
end

thresh = .05 / (size(alternate_ID_pct_bysub_byregion,2) / 2 *(size(alternate_ID_pct_bysub_byregion,2)-1));

[alternate_ID_region_correlations, alternate_ID_region_significance] = corrcoef(alternate_ID_pct_bysub_byregion_resid);

alternate_ID_region_correlations_thresh = alternate_ID_region_correlations .* (alternate_ID_region_significance < thresh);

figure
parcel_correlmat_figmaker(alternate_ID_region_correlations_thresh,assignments,[-.6 .6],'Residual Variable Region Correlations')
export_fig(gcf,'Residual Variable Region Correlations.pdf')

cifti_write_wHDR(alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta,[],'Variable_Regions_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta')
cifti_write_wHDR(alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval,[],'Variable_Regions_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval')
cifti_write_wHDR(alternate_ID_pct_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq,[],'Variable_Regions_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq')   

%% Variable region connectivity t-tests


percentregion_thresh = .5;

%all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8_selected.dtseries.nii');
variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');
all_variable_regions = sum(variable_regions,2);
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
mostcommon = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');

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

disp(alternatecount)
disp(alternatecount / length(subjects))

%variable_regions = variable_regions(:,1:end-1);

mainout = zeros(size(mostcommon,1),size(variable_regions,2));
altout = zeros(size(mostcommon,1),size(variable_regions,2));
tout = zeros(size(mostcommon,1),size(variable_regions,2));
cluscorrectedtout = zeros(size(mostcommon,1),size(variable_regions,2));
correctedtthresh = tinv(1-(.05 / ncortverts / size(variable_regions,2)),length(subjects)-2);
uncorrectedtthresh = tinv(1-(.001/ size(variable_regions,2)),length(subjects)-2);
for r = 1:size(variable_regions,2)
    
    %     alternate(:,r) = alternate(:,r) ./ alternatecount(r);
    %     main(:,r) = main(:,r) ./ (length(subjects) - alternatecount(r));
    
    mainout(1:ncortverts,r) = mean(main{r},2);
    altout(1:ncortverts,r) = mean(alternate{r},2);
    
    [H,P,CI,STATS] = ttest2(main{r}',alternate{r}');
    tout(1:ncortverts,r) = STATS.tstat;
    
    
    cluscorrectedtout(1:ncortverts,r) = Cluster_correct_ttest2_cifti(main{r},alternate{r},uncorrectedtthresh,.05,[]);
    
end



% cifti_write_wHDR(mainout,[],'MeanConnectivity_mainID_byregion')
% cifti_write_wHDR(altout,[],'MeanConnectivity_alternateID_byregion')
% cifti_write_wHDR(tout,[],['MeanConnectivity_T_MainVsAlternateID_byregion_correctedT' num2str(tthresh)])
% cifti_write_wHDR((mainout - altout),[],'MeanConnectivity_Diff_MainVsAlternateID_byregion')

cifti_write_wHDR(mainout,[],'MeanConnectivity_mainID_bywholeregion')
cifti_write_wHDR(altout,[],'MeanConnectivity_alternateID_bywholeregion')
cifti_write_wHDR(tout,[],['MeanConnectivity_T_MainVsAlternateID_bywholeregion_correctedT' num2str(correctedtthresh)])
cifti_write_wHDR((mainout - altout),[],'MeanConnectivity_Diff_MainVsAlternateID_bywholeregion')
cifti_write_wHDR(cluscorrectedtout,[],'MeanConnectivity_T_MainVsAlternateID_bywholeregion_cluscorrected')


%% Variable region network connectivity distributions

percentregion_thresh = .1;

%all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8_selected.dtseries.nii');
variable_regions = ft_read_cifti_mod('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii'); variable_regions = variable_regions.data;
all_variable_regions = sum(variable_regions,2);
networkconnection_bysub = ft_read_cifti_mod('Templatematch_dice_bysubject_kden0.05.dtseries.nii'); networkconnection_bysub = networkconnection_bysub.data;
mostcommon = ft_read_cifti_mod('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii'); mostcommon = mostcommon.data;

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

ncortverts = 59412;

%main = zeros(ncortverts,size(variable_regions,2));
%alternate = zeros(ncortverts,size(variable_regions,2));
alternatecount = zeros(1,size(variable_regions,2));

for s = 1:length(subjects)
    disp(['Subject ' num2str(s) ': ' subjects{s}])
    tmask = load(tmasks{s});
    subdata = ft_read_cifti_mod(ciftifiles{s});
    subdata = subdata.data(:,logical(tmask));
    
    for r = 1:size(variable_regions,2)
        
        if s==1
            alternate{r} = zeros(ncortverts,0);
            main{r} = zeros(ncortverts,0);
            
            alttomain{r} = [];
            alttoalt{r} = [];
            maintoalt{r} = [];
            maintomain{r} = [];
            
        end
        
        mainIDs_in_region = mostcommon(logical(variable_regions(:,r)));
        mainIDs_in_region(mainIDs_in_region<1) = [];
        mainID = mode(mainIDs_in_region); 
        altID = mean(all_variable_regions(logical(variable_regions(:,r))));
        
        altnetworkindices = find((networkconnection_bysub(:,s)==altID) & (~logical(variable_regions(:,r))));
        mainnetworkindices = find((networkconnection_bysub(:,s)==mainID) & (~logical(variable_regions(:,r))));
        
        if (nnz(networkconnection_bysub(logical(variable_regions(:,r)),s)==altID) / nnz(variable_regions(:,r))) >= percentregion_thresh
            
            altindices = find(logical(variable_regions(:,r)) & (networkconnection_bysub(:,s)==altID));
            
            
            
            alttomain{r}(end+1,1) = paircorr_mod(mean(subdata(altindices,:),1)',mean(subdata(mainnetworkindices,:),1)');
            alttoalt{r}(end+1,1) = paircorr_mod(mean(subdata(altindices,:),1)',mean(subdata(altnetworkindices,:),1)');
         end
         
         if (nnz(networkconnection_bysub(logical(variable_regions(:,r)),s)==mainID) / nnz(variable_regions(:,r))) >= percentregion_thresh
%        else
            
            mainindices = find(logical(variable_regions(:,r)) & (networkconnection_bysub(:,s)==mainID));
            
            maintomain{r}(end+1,1) = paircorr_mod(mean(subdata(mainindices,:),1)',mean(subdata(mainnetworkindices,:),1)');
            maintoalt{r}(end+1,1) = paircorr_mod(mean(subdata(mainindices,:),1)',mean(subdata(altnetworkindices,:),1)');
        end
            
         
    end
end

for r = 1:size(variable_regions,2)
    vec1 = FisherTransform(maintoalt{r});  vec1(isnan(vec1)) = [];
    vec2 = FisherTransform(alttoalt{r}); vec2(isnan(vec2)) = [];
    vecmin = min([vec1;vec2]); vecmax = max([vec1;vec2]);
    skew = skewness([vec1;vec2]);
    kurt = kurtosis([vec1;vec2]);
    nvals = length([vec1;vec2]);
    bimod_coeffs(r,1) = (skew^2 + 1) / (kurt + (3*(nvals-1)^2) / ((nvals-2)*(nvals-3)));
    %[dip, hard_p_value(r,1), xlow,xup]=HartigansDipSignifTest(sort([vec1;vec2]'),1000);
    
    vec1(end+1:length(vec2)) = NaN; vec2(end+1:length(vec1)) = NaN;
    figure;
    [N,C] = hist([vec1 vec2],15);
    bar(C,N,2)
    title(gca,['Region ' num2str(r) ', to Alternate System'])
    legend(gca,'Main','Alt')
    set(gcf,'Color','white')
    set(gca,'FontSize',20)
end


    
 for r = 1:size(variable_regions,2)
    
    vec1 = FisherTransform(maintomain{r});  vec1(isnan(vec1)) = [];
    vec2 = FisherTransform(alttomain{r});  vec2(isnan(vec2)) = [];
    vecmin = min([vec1;vec2]); vecmax = max([vec1;vec2]);
    skew = skewness([vec1;vec2]);
    kurt = kurtosis([vec1;vec2]);
    bimod_coeffs(r,2) = (skew^2 + 1) / (kurt + (3*(nvals-1)^2) / ((nvals-2)*(nvals-3)));
    %[dip, hard_p_value(r,2), xlow,xup]=HartigansDipSignifTest(sort([vec1;vec2]'),1000);
    
    vec1(end+1:length(vec2)) = NaN; vec2(end+1:length(vec1)) = NaN;
    figure;
    [N,C] = hist([vec1 vec2],15);
    bar(C,N,2)
    title(gca,['Region ' num2str(r) ', to Main System'])
    legend(gca,'Main','Alt')
    set(gcf,'Color','white')
    set(gca,'FontSize',20)
end

%% Variable region connectivity t-tests, confound regressed


percentregion_thresh = .5;

%all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8_selected.dtseries.nii');
variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');
all_variable_regions = sum(variable_regions,2);
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
mostcommon = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

ncortverts = 59412;

%main = zeros(ncortverts,size(variable_regions,2));
%alternate = zeros(ncortverts,size(variable_regions,2));
alternatecount = zeros(1,size(variable_regions,2));

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
sulcdepth_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_L = sulcdepth_group_L.cdata(logical(maskL));
curv_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.curvature.32k_fs_LR.shape.gii'); curv_group_L = curv_group_L.cdata(logical(maskL));

maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
sulcdepth_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_R = sulcdepth_group_R.cdata(logical(maskR));
curv_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.curvature.32k_fs_LR.shape.gii'); curv_group_R = curv_group_R.cdata(logical(maskR));

sulc_within_region = zeros(length(subjects),size(variable_regions,2));
curv_within_region = zeros(length(subjects),size(variable_regions,2));

connectivitymaps = zeros(ncortverts,length(subjects),size(variable_regions,2));

alternate_binary = zeros(length(subjects),size(variable_regions,2));

for s = 1:length(subjects)
    disp(subjects{s})
    tmask = load(tmasks{s});
    subdata = ft_read_cifti_mod(ciftifiles{s});
    subdata = subdata.data(:,logical(tmask));
    
    sulcdepth_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_diff = [sulcdepth_L.cdata(logical(maskL)) - sulcdepth_group_L ; sulcdepth_R.cdata(logical(maskR)) - sulcdepth_group_R];
    
    curv_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.curvature.32k_fs_LR.shape.gii']);
    curv_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.curvature.32k_fs_LR.shape.gii']);
    curv_diff = [curv_L.cdata(logical(maskL)) - curv_group_L ; curv_R.cdata(logical(maskR)) - curv_group_R];
    
    
    AD_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.ArealDistortion_32k_fs_LR.shape.gii']);
    AD_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.ArealDistortion_32k_fs_LR.shape.gii']);
    AD = abs([AD_L.cdata(logical(maskL)) ; AD_R.cdata(logical(maskR))]);
    
    subSNR = ft_read_cifti_mod(['/data/cn4/evan/RestingState/Ind_variability/120_108_SNR/' subjects{s} '_SNR.dtseries.nii']);
    SNR = subSNR.data(1:ncortverts,:);
    
    
    tmask = load(tmasks{s});
    ntimepoints(s,1) = nnz(tmask);
    
    fd_file = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/FCPROCESS_bandpass_interp_nosmooth/' subjects{s} '/total_FD.txt'];
    if exist(fd_file)
        FD = load(fd_file);
        meanFD(s,1) = mean(FD(logical(tmask)));
    else
        fd_file = ['/data/cn5/selfRegulation/V4Process/initial_process/' subjects{s} '/' subjects{s} '/total_FD.txt'];
        FD = load(fd_file);
        meanFD(s,1) = mean(FD(logical(tmask)));
    end
    
    for r = 1:size(variable_regions,2)
        
        sulc_within_region(s,r) = mean(abs(sulcdepth_diff(logical(variable_regions(:,r)))));
        curv_within_region(s,r) = mean(abs(curv_diff(logical(variable_regions(:,r)))));
        AD_within_region(s,r) = mean(AD(logical(variable_regions(:,r))));
        SNR_within_region(s,r) = mean(SNR(logical(variable_regions(:,r))));
        
        %if s==1
            %alternate{r} = zeros(ncortverts,0);
            %main{r} = zeros(ncortverts,0);
            
        %end
        
        mainID = mode(mostcommon(logical(variable_regions(:,r)))); 
        altID = mean(all_variable_regions(logical(variable_regions(:,r))));
        if (nnz(networkconnection_bysub(logical(variable_regions(:,r)),s)==altID) / nnz(variable_regions(:,r))) >= percentregion_thresh
            alternatecount(r) = alternatecount(r)+1;
            alternate_binary(s,r) = 1;
        end
            
            indices = logical(variable_regions(:,r));% .* (networkconnection_bysub(:,s)==altID));
            
            corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
            corrpattern(isnan(corrpattern)) = 0;
            
            connectivitymaps(:,s,r) = corrpattern(1:ncortverts);
            
            %alternate(:,r) = alternate(:,r) + corrpattern(1:ncortverts);
            %alternate{r}(:,end+1) = corrpattern(1:ncortverts);
            
%         else
%             
%             indices = logical(variable_regions(:,r));% .* (networkconnection_bysub(:,s)==mainID));
%             
%             corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
%             corrpattern(isnan(corrpattern)) = 0;
%             
%             %main(:,r) = main(:,r) + corrpattern(1:ncortverts);
%             %main{r}(:,end+1) = corrpattern(1:ncortverts);
%         end
    end
end

disp(alternatecount)
disp(alternatecount / length(subjects))

%variable_regions = variable_regions(:,1:end-1);

mainout = zeros(size(mostcommon,1),size(variable_regions,2));
altout = zeros(size(mostcommon,1),size(variable_regions,2));
tout = zeros(size(mostcommon,1),size(variable_regions,2));
cluscorrectedtout = zeros(size(mostcommon,1),size(variable_regions,2));
correctedtthresh = tinv(1-(.05 / ncortverts / size(variable_regions,2)),length(subjects)-2);
uncorrectedtthresh = tinv(1-(.001/ size(variable_regions,2)),length(subjects)-2);
for r = 1:size(variable_regions,2)
    
    connectivitymaps_regressed = zeros(size(mostcommon,1),length(subjects));
    
    for vert = 1:size(connectivitymaps,1)
        STATS = regstats(squeeze(connectivitymaps(vert,:,r))',[sulc_within_region(:,r) curv_within_region(:,r) ntimepoints meanFD AD_within_region(:,r) SNR_within_region(:,r)],'linear',{'r'});
        connectivitymaps_regressed(vert,:) = STATS.r';
    end
    
    %     alternate(:,r) = alternate(:,r) ./ alternatecount(r);
    %     main(:,r) = main(:,r) ./ (length(subjects) - alternatecount(r));
    
    mainout(1:ncortverts,r) = mean(connectivitymaps_regressed(1:ncortverts,logical(alternate_binary(:,r)==0)),2);
    altout(1:ncortverts,r) =  mean(connectivitymaps_regressed(1:ncortverts,logical(alternate_binary(:,r)==1)),2);
    
    [H,P,CI,STATS] = ttest2(connectivitymaps_regressed(:,logical(alternate_binary(:,r)==0))',connectivitymaps_regressed(:,logical(alternate_binary(:,r)==1))');
    tout(:,r) = STATS.tstat;
    
    
    cluscorrectedtout(:,r) = Cluster_correct_ttest2_cifti(connectivitymaps_regressed(:,logical(alternate_binary(:,r)==0)),connectivitymaps_regressed(:,logical(alternate_binary(:,r)==1)),uncorrectedtthresh,.05,[]);
    
end



% cifti_write_wHDR(mainout,[],'MeanConnectivity_mainID_byregion')
% cifti_write_wHDR(altout,[],'MeanConnectivity_alternateID_byregion')
% cifti_write_wHDR(tout,[],['MeanConnectivity_T_MainVsAlternateID_byregion_correctedT' num2str(tthresh)])
% cifti_write_wHDR((mainout - altout),[],'MeanConnectivity_Diff_MainVsAlternateID_byregion')

cifti_write_wHDR(mainout,[],'MeanConnectivityNewregressed_mainID_bywholeregion')
cifti_write_wHDR(altout,[],'MeanConnectivityNewregressed_alternateID_bywholeregion')
cifti_write_wHDR(tout,[],['MeanConnectivityNewregressed_T_MainVsAlternateID_bywholeregion_correctedT' num2str(correctedtthresh)])
cifti_write_wHDR((mainout - altout),[],'MeanConnectivityNewregressed_Diff_MainVsAlternateID_bywholeregion')
cifti_write_wHDR(cluscorrectedtout,[],'MeanConnectivityNewregressed_T_MainVsAlternateID_bywholeregion_cluscorrected')

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

%% Cluster Movement Metrics

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');
subjects = subjects(121:end);%subjects(1:120);

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
groupavgsystems = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');%'Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');

orig_avgsystems = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR.dtseries.nii');
networkconnection_bysub(logical((orig_avgsystems==13) + (orig_avgsystems==14)),:) = 0;

mask{1} = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'); mask{1} = ~mask{1}.cdata;
mask{2} = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'); mask{2} = ~mask{2}.cdata;

% verts{1} = gifti('/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii'); verts{1} = verts{1}.vertices;
% verts{2} = gifti('/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii'); verts{2} = verts{2}.vertices;


nsurfverts = nnz(mask{1}) + nnz(mask{2});

neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');
load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
distances = distances(1:nsurfverts,1:nsurfverts);

pct_overlap = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
size_change = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
centroid_movement = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
border_movement = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;
border_movement_directed = ones(nsurfverts,size(networkconnection_bysub,2)) * NaN;



groupclusters_successfully_matched = zeros(nsurfverts,1);

for s = 1:size(networkconnection_bysub,2)
    disp(['Subject ' num2str(s)])
    fprintf('   ')
    [groupmatches, submatches] = networkClusters_matchtogroup(networkconnection_bysub(:,s),distances,groupavgsystems);
    
    groupclusters_successfully_matched = groupclusters_successfully_matched + (sum(logical(groupmatches),2) / size(networkconnection_bysub,2));
    
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
        this_size_change = (sum(SA_LR(logical(submatches(:,i)))) - sum(SA_LR(logical(groupmatches(:,i))))) ./ sum(SA_LR(logical(groupmatches(:,i))));
        
        
        size_change(logical(groupmatches(:,i)),s) = this_size_change;
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

out = zeros(size(groupavgsystems));
out(1:nsurfverts) = groupclusters_successfully_matched;
cifti_write_wHDR(out,[],'Group_Cluster_Match_Pct')

out = zeros(size(groupavgsystems));
out(1:nsurfverts) = nanmean(border_movement,2);
cifti_write_wHDR(out,[],'Avg_Cluster_Border_Movement')

out = zeros(size(groupavgsystems));
out(1:nsurfverts) = nanmean(abs(size_change),2);
cifti_write_wHDR(out,[],'Avg_Cluster_Size_Change')

out = zeros(size(groupavgsystems,1),size(networkconnection_bysub,2));
out(1:nsurfverts,:) = size_change;
cifti_write_wHDR(out,[],'Cluster_Size_Change_bysub')

SAs = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned_clusterSAs.dtseries.nii');
std_size = nanstd((out .* repmat(SAs,1,size(out,2)) + repmat(SAs,1,size(out,2))) ./ repmat(SAs,1,size(out,2)),[],2);
cifti_write_wHDR(std_size,[],'Std_Cluster_Sizes.dtseries.nii')

out = zeros(size(groupavgsystems,1),size(networkconnection_bysub,2));
out(1:nsurfverts,:) = border_movement_directed;
cifti_write_wHDR(out,[],'Border_Movement_bysub')

%save('border_distance_directed_bysub.mat','border_movement_directed','-v7.3')



%% Border Segment Movement Correlations

% ncortverts = 59412;
% 
% border_movement_directed = cifti_read('Border_Movement_bysub.dtseries.nii');
% border_movement_directed = border_movement_directed(1:ncortverts,:);
% 
% 
% neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');
% 
% mostcommon_name = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii';
% 
% mostcommon = cifti_read(mostcommon_name);
% 
% 
% edges = zeros(size(mostcommon));
% 
% for i = 1:ncortverts
%     
%     vertneighs = neighbors(i,2:end); vertneighs(isnan(vertneighs)) = [];
%     if any(mostcommon(vertneighs) == mostcommon(i))
%         edges(i) = mostcommon(i);
%     end
% end
% 
% %cifti_write_wHDR(edges,[],'Variability_LR_consensus_dice_mostcommon_clean50mm_edges')
% 
% IDs = unique(mostcommon); IDs(IDs==0) = [];
% 
% segments = zeros(size(edges,1),0);
% 
% segmentIDs = zeros(0,2);
% 
% segment_movement = zeros(0,size(border_movement_directed,2));
% 
% for ID1 = IDs(:)'
%     for ID2 = IDs(:)'
%         if ID1 ~= ID2
%             
%             ID1_nexttoID2_edges = zeros(size(edges));
%             
%             ID1_edgeverts = find(edges==ID1);
%             
%             for vert = ID1_edgeverts(:)'
%                 
%                 vertneighs = neighbors(vert,2:end); vertneighs(isnan(vertneighs)) = [];
%                 
%                 if any(edges(vertneighs)==ID2)
%                     
%                     ID1_nexttoID2_edges(vert) = 1;
%                 end
%             end
%             
%             temp = metric_cluster_cifti(ID1_nexttoID2_edges,.5,1.5,0);
%             
%             segments(:,end+1:end+size(temp,2)) = temp;
%             
%             segmentIDs(end+1:end+size(temp,2),:) = repmat([ID1 ID2],size(temp,2),1);
%             
%             for segnum = 1:size(temp,2)
%                 segment_movement(end+1,:) = mean(border_movement_directed(logical(temp(:,segnum)),:),1);
%             end
%             
%         end
%     end
% end
% 
% segment_movement_correlations = zeros(size(segment_movement,1));
% segment_movement_significance = zeros(size(segment_movement,1));
% for i = 1:size(segment_movement,1)
%     for j = 1:size(segment_movement,1)
%         inds = logical((~isnan(segment_movement(i,:))) .* (~isnan(segment_movement(j,:))));
%         
%         [R, P] = corrcoef(segment_movement(i,inds)',segment_movement(j,inds)');
%         
%         segment_movement_correlations(i,j) = R(1,2);
%         segment_movement_significance(i,j) = P(1,2);
%     end
% end
% %thresh = .05;
% thresh = .05 / (length(segment_movement_correlations) / 2 *(length(segment_movement_correlations)-1));
% segment_movement_correlations = segment_movement_correlations .* (segment_movement_significance < thresh);
% 
% figure
% parcel_correlmat_figmaker(segment_movement_correlations,segmentIDs(:,1),[-.6 .6],'Border Movement Correlations')
% export_fig(gcf,'Border Movement Correlations by Segment.pdf')



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


%% Cluster Size Correlations
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
groupavgsystems = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');%'Variability_LR_consensus_dice_mostcommon_clean50mm.dtseries.nii');
size_change = cifti_read('Cluster_Size_Change_bysub.dtseries.nii');
%cluster_SAs = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned_clusterSAs.dtseries.nii');
cluster_SA_thresh = 100;

IDs = unique(groupavgsystems); IDs(IDs==0) = [];
groupavgsystems_clusters = zeros(length(groupavgsystems),0);
for ID = IDs(:)'
    outputcifti = metric_cluster_cifti_surfacearea(groupavgsystems,ID-.01,ID+.01,cluster_SA_thresh);
    groupavgsystems_clusters(:,end+1:end+size(outputcifti,2)) = outputcifti;
end


%pct_overlap_byclus = zeros(size(networkconnection_bysub,2),size(groupavgsystems_clusters,2));
%centroid_movement_byclus = zeros(size(networkconnection_bysub,2),size(groupavgsystems_clusters,2));
%border_movement_byclus = zeros(size(networkconnection_bysub,2),size(groupavgsystems_clusters,2));
size_change_byclus = zeros(size(networkconnection_bysub,2),size(groupavgsystems_clusters,2));
for i = 1:size(groupavgsystems_clusters,2)
    %pct_overlap_byclus(:,i) = mean(pct_overlap(logical(groupavgsystems_clusters(:,i)),:),1);
    %centroid_movement_byclus(:,i) = mean(centroid_movement(logical(groupavgsystems_clusters(:,i)),:),1);
    %border_movement_byclus(:,i) = nanmean(border_movement(logical(groupavgsystems_clusters(:,i)),:),1);
    size_change_byclus(:,i) = nanmean(size_change(logical(groupavgsystems_clusters(:,i)),:),1);
end

size_change_rs = zeros(size(groupavgsystems_clusters,2),size(groupavgsystems_clusters,2));
size_change_ps = size_change_rs;
%pct_overlap_rs = size_change_rs; pct_overlap_ps = pct_overlap_rs; centroid_movement_rs = pct_overlap_rs; centroid_movement_ps = pct_overlap_rs; border_movement_rs = pct_overlap_rs; border_movement_ps = pct_overlap_rs; size_change_rs = pct_overlap_rs; size_change_ps = pct_overlap_rs;

count_matched = zeros(size(groupavgsystems));

for i = 1:size(groupavgsystems_clusters,2)
    nonnaninds_i = ~isnan(size_change_byclus(:,i));
    count_matched(logical(groupavgsystems_clusters(:,i))) = nnz(nonnaninds_i);
    for j = 1:size(groupavgsystems_clusters,2)
        if i ~= j
            nonnaninds_j = ~isnan(size_change_byclus(:,j));
            nonnaninds_ij = logical(nonnaninds_i .* nonnaninds_j);
            if nnz(nonnaninds_ij) > 2
                %[pct_overlap_rs(i,j), pct_overlap_ps(i,j)] = corr(pct_overlap_byclus(nonnaninds_ij,i),pct_overlap_byclus(nonnaninds_ij,j));
                %[centroid_movement_rs(i,j), centroid_movement_ps(i,j)] = corr(centroid_movement_byclus(nonnaninds_ij,i),centroid_movement_byclus(nonnaninds_ij,j));
                %[border_movement_rs(i,j), border_movement_ps(i,j)] = corr(border_movement_byclus(nonnaninds_ij,i),border_movement_byclus(nonnaninds_ij,j));
                [size_change_rs(i,j), size_change_ps(i,j)] = corr(size_change_byclus(nonnaninds_ij,i),size_change_byclus(nonnaninds_ij,j));
            end
        end
    end
end

thresh = .05 / nnz(tril(ones(size(size_change_rs)),1));

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

assignments = zeros(size(groupavgsystems_clusters,2),1);
for i = 1:length(assignments)
    assignments(i) = mean(groupavgsystems(logical(groupavgsystems_clusters(:,i))));
end

% figure
% parcel_correlmat_figmaker((pct_overlap_rs .* (pct_overlap_ps < thresh)),assignments,[-.6 .6],'Percent Overlap Correlations')
% export_fig(gcf,'Percent Overlap Correlations.pdf')
% 
% figure
% parcel_correlmat_figmaker((centroid_movement_rs .* (centroid_movement_ps < thresh)),assignments,[-.6 .6],'Centroid Movement Correlations')
% export_fig(gcf,'Centroid Movement Correlations.pdf')
% 
% figure
% parcel_correlmat_figmaker((border_movement_rs .* (border_movement_ps < thresh)),assignments,[-.6 .6],'Border Movement Correlations')
% export_fig(gcf,'Border Movement Correlations.pdf')


parcel_correlmat_figmaker((size_change_rs .* (size_change_ps < thresh)),assignments,[-.6 .6],'Size Change Correlations')
export_fig(gcf,'Size Change Correlations.pdf')

  
%% Cluster Size Regressing data quality

ncortverts = 59412;

mostcommon_name = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii';
toosmall_clusnums = [9 17 19 24 43 46 52 59 75 79];

mostcommon = cifti_read(mostcommon_name);
IDs = unique(mostcommon); IDs(IDs==0) = [];
mostcommon_clusters = zeros(length(mostcommon),0);
for ID = IDs(:)'
    outputcifti = metric_cluster_cifti_surfacearea(mostcommon,ID-.01,ID+.01,10);
    mostcommon_clusters(:,end+1:end+size(outputcifti,2)) = outputcifti;
end
mostcommon_clusters_use = ~any(mostcommon_clusters(:,toosmall_clusnums),2);

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
AD = zeros(ncortverts,length(subjects));
SNR = zeros(ncortverts,length(subjects));
AD_clus = zeros(size(mostcommon_clusters,2),length(subjects));
SNR_clus = zeros(size(mostcommon_clusters,2),length(subjects));

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
sulcdepth_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_L = sulcdepth_group_L.cdata(logical(maskL));
curv_group_L = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.curvature.32k_fs_LR.shape.gii'); curv_group_L = curv_group_L.cdata(logical(maskL));

maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
sulcdepth_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.sulc.32k_fs_LR.shape.gii'); sulcdepth_group_R = sulcdepth_group_R.cdata(logical(maskR));
curv_group_R = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.curvature.32k_fs_LR.shape.gii'); curv_group_R = curv_group_R.cdata(logical(maskR));


for s = 1:length(subjects)
    disp(num2str(s))
    sulcdepth_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.sulc.32k_fs_LR.shape.gii']);
    sulcdepth_diff(:,s) = [sulcdepth_L.cdata(logical(maskL)) - sulcdepth_group_L ; sulcdepth_R.cdata(logical(maskR)) - sulcdepth_group_R];
        
    curv_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.curvature.32k_fs_LR.shape.gii']);
    curv_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.curvature.32k_fs_LR.shape.gii']);
    curv_diff(:,s) = [curv_L.cdata(logical(maskL)) - curv_group_L ; curv_R.cdata(logical(maskR)) - curv_group_R];
    
    AD_L = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.L.ArealDistortion_32k_fs_LR.shape.gii']);
    AD_R = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.R.ArealDistortion_32k_fs_LR.shape.gii']);
    AD(:,s) = abs([AD_L.cdata(logical(maskL)) ; AD_R.cdata(logical(maskR))]);
    
    subSNR = ft_read_cifti_mod(['/data/cn4/evan/RestingState/Ind_variability/120_108_SNR/' subjects{s} '_SNR.dtseries.nii']);
    SNR(:,s) = subSNR.data(1:59412,:);
    
    
    for clusnum = 1:size(mostcommon_clusters,2)
        sulcdepth_diff_clus(clusnum,s) = mean(abs(sulcdepth_diff(logical(mostcommon_clusters(:,clusnum)),s)));
        curv_diff_clus(clusnum,s) = mean(abs(curv_diff(logical(mostcommon_clusters(:,clusnum)),s)));
        AD_clus(clusnum,s) = mean(abs(AD(logical(mostcommon_clusters(:,clusnum)),s)));
        SNR_clus(clusnum,s) = mean(abs(SNR(logical(mostcommon_clusters(:,clusnum)),s)));
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

size_change = cifti_read('Cluster_Size_Change_bysub.dtseries.nii');

Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta = zeros(66697,6);
Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval = zeros(66697,6);
Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq = zeros(66697,1);
Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_resid = ones(66697,size(size_change,2)) * NaN;

verts = find(any(~isnan(size_change),2) .* any(size_change > 0 , 2));
for vert = verts(:)'
    
    clusnum_ind = find(mostcommon_clusters(vert,:));
    
    inds = logical(~isnan(size_change(vert,:)));
    
    STATS = regstats(size_change(vert,inds)',[sulcdepth_diff_clus(clusnum_ind,inds)' curv_diff_clus(clusnum_ind,inds)' ntimepoints(inds) meanFD(inds) AD_clus(clusnum_ind,inds)' SNR_clus(clusnum_ind,inds)'],'linear',{'tstat','r','adjrsquare'});
    Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta(vert,:) = STATS.tstat.beta(2:end);
    Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval(vert,:) = STATS.tstat.pval(2:end);
    Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq(vert) = STATS.adjrsquare;
    Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_resid(vert,inds) = STATS.r;
    
end

cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta.*repmat(mostcommon_clusters_use,1,6),[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_beta')
cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval.*repmat(mostcommon_clusters_use,1,6),[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval')
cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval * (size(mostcommon_clusters,2) - length(toosmall_clusnums)) .*repmat(mostcommon_clusters_use,1,6),[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_pval_corrected')
cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq.*mostcommon_clusters_use,[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_adjrsq')
cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_resid,[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_resid')
%cifti_write_wHDR(Size_V_Sulc_Curv_Ntimepoints_FD_resid .* sign(cifti_read('Cluster_Size_Change_bysub.dtseries.nii')),[],'Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_resid_directed')



% border_shift = abs(cifti_read('Border_Movement_bysub.dtseries.nii'));
% 
% Border_V_Sulc_Curv_Ntimepoints_FD_beta = zeros(66697,4);
% Border_V_Sulc_Curv_Ntimepoints_FD_pval = zeros(66697,4);
% Border_V_Sulc_Curv_Ntimepoints_FD_resid = ones(66697,size(border_shift,2)) * NaN;
% 
% verts = find(any(~isnan(border_shift),2) .* any(border_shift > 0 , 2));
% 
% for vert = verts(:)'
%     
%     inds = logical(~isnan(border_shift(vert,:)));
%     
%     STATS = regstats(border_shift(vert,inds)',[sulcdepth_diff(vert,inds)' curv_diff(vert,inds)' ntimepoints(inds) meanFD(inds)],'linear',{'tstat','r'});
%     Border_V_Sulc_Curv_Ntimepoints_FD_beta(vert,:) = STATS.tstat.beta(2:end);
%     Border_V_Sulc_Curv_Ntimepoints_FD_pval(vert,:) = STATS.tstat.pval(2:end);
%     Border_V_Sulc_Curv_Ntimepoints_FD_resid(vert,inds) = STATS.r;
%     
% end
% 
% cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_beta,[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_beta')
% cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_pval,[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_pval')
% cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_resid,[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid')
% cifti_write_wHDR(Border_V_Sulc_Curv_Ntimepoints_FD_resid .* sign(cifti_read('Border_Movement_bysub.dtseries.nii')),[],'BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid_directed')



%% Residual Cluster Size and Border Movement correlations

ncortverts = 59412;
cluster_SA_thresh = 100;

mostcommon_name = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii';
mostcommon = cifti_read(mostcommon_name);
IDs = unique(mostcommon); IDs(IDs==0) = [];
mostcommon_clusters = zeros(length(mostcommon),0);
for ID = IDs(:)'
outputcifti = metric_cluster_cifti_surfacearea(mostcommon,ID-.01,ID+.01,cluster_SA_thresh);
mostcommon_clusters(:,end+1:end+size(outputcifti,2)) = outputcifti;
end

border_movement_resid = cifti_read('BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid_directed.dtseries.nii');
%size_change_resid = cifti_read('Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_resid_directed.dtseries.nii');
size_change_resid = cifti_read('Cluster_Size_V_Sulc_Curv_Ntimepoints_FD_AD_SNR_resid.dtseries.nii');

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

% figure
% parcel_correlmat_figmaker((border_movement_resid_rs .* (border_movement_resid_ps < thresh)),assignments,[-.6 .6],'Residual Border Movement Correlations')
% export_fig(gcf,'Residual Border Movement Correlations.pdf')

figure
parcel_correlmat_figmaker((size_change_resid_rs .* (size_change_resid_ps < thresh)),assignments,[-.6 .6],'Residual Size Change Correlations')
export_fig(gcf,'Residual Size Change Correlations.pdf')

%% Residual Border Segment movement correlations

% ncortverts = 59412;
% 
% border_movement_directed = cifti_read('BorderMovement_V_Sulc_Curv_Ntimepoints_FD_resid_directed.dtseries.nii');
% border_movement_directed = border_movement_directed(1:ncortverts,:);
% 
% 
% neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');
% 
% mostcommon_name = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii';
% 
% mostcommon = cifti_read(mostcommon_name);
% 
% edges = zeros(size(mostcommon));
% 
% for i = 1:ncortverts
%     
%     vertneighs = neighbors(i,2:end); vertneighs(isnan(vertneighs)) = [];
%     if any(mostcommon(vertneighs) == mostcommon(i))
%         edges(i) = mostcommon(i);
%     end
% end
% 
% %cifti_write_wHDR(edges,[],'Variability_LR_consensus_dice_mostcommon_clean50mm_edges')
% 
% IDs = unique(mostcommon); IDs(IDs==0) = [];
% 
% segments = zeros(size(edges,1),0);
% 
% segmentIDs = zeros(0,2);
% 
% segment_movement = zeros(0,size(border_movement_directed,2));
% 
% for ID1 = IDs(:)'
%     for ID2 = IDs(:)'
%         if ID1 ~= ID2
%             
%             ID1_nexttoID2_edges = zeros(size(edges));
%             
%             ID1_edgeverts = find(edges==ID1);
%             
%             for vert = ID1_edgeverts(:)'
%                 
%                 vertneighs = neighbors(vert,2:end); vertneighs(isnan(vertneighs)) = [];
%                 
%                 if any(edges(vertneighs)==ID2)
%                     
%                     ID1_nexttoID2_edges(vert) = 1;
%                 end
%             end
%             
%             temp = metric_cluster_cifti(ID1_nexttoID2_edges,.5,1.5,0);
%             
%             segments(:,end+1:end+size(temp,2)) = temp;
%             
%             segmentIDs(end+1:end+size(temp,2),:) = repmat([ID1 ID2],size(temp,2),1);
%             
%             for segnum = 1:size(temp,2)
%                 segment_movement(end+1,:) = mean(border_movement_directed(logical(temp(:,segnum)),:),1);
%             end
%             
%         end
%     end
% end
% 
% segment_movement_correlations = zeros(size(segment_movement,1));
% segment_movement_significance = zeros(size(segment_movement,1));
% for i = 1:size(segment_movement,1)
%     for j = 1:size(segment_movement,1)
%         inds = logical((~isnan(segment_movement(i,:))) .* (~isnan(segment_movement(j,:))));
%         
%         [R, P] = corrcoef(segment_movement(i,inds)',segment_movement(j,inds)');
%         
%         segment_movement_correlations(i,j) = R(1,2);
%         segment_movement_significance(i,j) = P(1,2);
%     end
% end
% %thresh = .05;
% thresh = .05 / (length(segment_movement_correlations) / 2 *(length(segment_movement_correlations)-1));
% segment_movement_correlations = segment_movement_correlations .* (segment_movement_significance < thresh);
% 
% figure
% parcel_correlmat_figmaker(segment_movement_correlations,segmentIDs(:,1),[-.6 .6],'Residual Border Movement Correlations')
% export_fig(gcf,'Residual Border Movement Correlations by Segment.pdf')




