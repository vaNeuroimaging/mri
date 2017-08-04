maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;%maskL = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/L.atlasroi_group_noproj.func.gii');maskL = maskL.cdata;%
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;%maskR = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/R.atlasroi_group_noproj.func.gii');maskR = maskR.cdata;%

surfcoordsL = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii'); surfcoordsL = surfcoordsL.vertices;
surfcoordsR = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii'); surfcoordsR = surfcoordsR.vertices;
surfcoords = [surfcoordsL ; surfcoordsR];

subcort_coords = load('/data/hcp-zfs/shared-nil/laumannt/120_parcellation/modified_cifti_network/cifti_normalwall/120/BOTH/cifti_coords_LR.txt');
subcort_coords = subcort_coords(59413:end,:);

T4_711_to_MNI_mat = [0.953903 -0.003872 -0.021291; -0.010402  0.950932 -0.016412; 0.020852  0.055067  0.946383];

for i = 1:size(subcort_coords,1)
    subcort_coords(i,:) = subcort_coords(i,:) * T4_711_to_MNI_mat;
end


distmat_use = zeros(length(maskL) + length(maskR) + length(subcort_coords));

distancesL = smartload('/data/cn4/evan/fsaverage_LR32k/Surface_distances_L.mat');
distancesR = smartload('/data/cn4/evan/fsaverage_LR32k/Surface_distances_R.mat');

distmat_use(1:length(maskL),1:length(maskL)) = distancesL;
distmat_use((1+length(maskL)):(length(maskL)+length(maskR)),(1+length(maskL)):(length(maskL)+length(maskR))) = distancesR;
distmat_use(1:length(maskL),(1+length(maskL)):(length(maskL)+length(maskR))) = 1000;%distancesL;
distmat_use((1+length(maskL)):(length(maskL)+length(maskR)),1:length(maskL)) = 1000;%distancesR;



distmat_use(1:(length(maskL)+length(maskR)),(length(maskL)+length(maskR)+1) : end) = pdist2(surfcoords,subcort_coords);
distmat_use((length(maskL)+length(maskR)+1) : end,1:(length(maskL)+length(maskR))) = pdist2(subcort_coords,surfcoords);
distmat_use((length(maskL)+length(maskR)+1) : end,(length(maskL)+length(maskR)+1) : end) = pdist2(subcort_coords,subcort_coords);

maskvec = logical([maskL;maskR;ones(length(subcort_coords),1)]);

distmat_use = distmat_use(maskvec,maskvec);

save('distmat_normalwall_surf_geodesic_vol_euc.mat','distmat_use','-v7.3')
%save('distmat_smallwall_homo_surf_geodesic_vol_euc.mat','distmat_use','-v7.3')

