load /data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_LR_infomap/corrmat.mat

assigns = load('/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_LR_infomap/watershed_Tk0005to005in0001_S1to1_surfxd20_INFMAP/120_subsurf_LR_minsize5_consensus.txt');
communities = unique(assigns(assigns>0));

parcels = cifti_read('/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_LR_nosmooth_watershedmerge_0.4_tweaked.dtseries.nii');
parcelIDs = unique(parcels(parcels>0)); 
threshs = [.01 : .002 : .03];

output = zeros(size(parcels,1),length(threshs));
ciftitemplatefile = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];

for threshnum = 1:length(threshs)
    thresh = threshs(threshnum);
    
    [threshmat r kden] = matrix_thresholder(all_water_corrmat,thresh,'kden');
    
    for parcelnum = 1:length(assigns);
        
        nodeconnections = logical(threshmat(:,parcelnum));
        degree = nnz(nodeconnections);
        if degree == 0
            participationcoeff = 0;
        else
            ratios = zeros(length(communities),1);
            
            for communitynum = 1:length(communities)
                numthiscommunityconnections = nnz(assigns(nodeconnections)==communities(communitynum));
                ratios(communitynum) = (numthiscommunityconnections/degree)^2;
            end
            
            participationcoeff = 1-sum(ratios);
        end
        output(parcels==parcelIDs(parcelnum),threshnum) = participationcoeff;
    end
end

output(parcels==0,:) = -1;
cifti_write_wHDR(output,ciftitemplatefile,'120_subsurf_LR_minsize5_consensus_participationcoeff')

