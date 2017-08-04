mapnumbers = [7 11 14 26 34 35 50 76 77 81 86 142 148 180 183 262];
parcelnumbers = {[1 26],75,237,251,197,39,257,[74 110],[236 273],180,21,[28 27 22 34 21 30],[7 74 110],264,99,154};

%maps = ft_read_cifti_mod('Cluster_probability_maps_sorted_10mm_40sub.dscalar.nii');
parcels = ft_read_cifti_mod('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii'); 

neighbors = cifti_neighbors('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii');

parcelsout = zeros(size(parcels.data,1),length(mapnumbers));

for i = 1:length(mapnumbers)
    for j = parcelnumbers{i}
        parcelinds = find(parcels.data==j);
        parcelneighs = unique(neighbors(parcelinds,2:7)); parcelneighs(isnan(parcelneighs)) = [];
        parcelborders = setdiff(parcelneighs,parcelinds);
        parcelsout(parcelborders,i) = 1;
    end
end

parcels.data = parcelsout;
ft_write_cifti_mod('Cluster_probability_maps_sorted_10mm_parcelsmatch',parcels)