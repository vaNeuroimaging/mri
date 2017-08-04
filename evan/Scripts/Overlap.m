parcels = cifti_read('../../120_subsurf_LR_nosmooth_watershedmerge_0.4_tweaked_gooddata.dtseries.nii');
parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];

Infomap = cifti_read('120_subsurf_LR_minsize5_consensus.dtseries.nii');

Power = cifti_read('/data/cn4/evan/ROIs/Power_Neuron_2011_LR.func.gii.dtseries.nii');

output = ones(66697,1) * -1;

for IDnum = 1:length(parcelIDs)
ID = parcelIDs(IDnum);
parcelinds = find(parcels==ID);
if any(Power(parcelinds)>0) && any(Infomap(parcelinds)>0)
    nonzeroparcelinds = intersect(intersect(parcelinds,find(Power>0)),find(Infomap>0));
overlap(IDnum) = nnz(Power(nonzeroparcelinds)==Infomap(nonzeroparcelinds(1))) / numel(nonzeroparcelinds);
output(parcelinds) = overlap(IDnum);
else
overlap(IDnum) = NaN;
output(parcelinds) = -1;
end
end

disp(['Overlap by parcel: ' num2str(nanmean(overlap))])


allparcelinds = intersect(find(Power>0),find(newInfomap>0));
areaoverlap = nnz(Infomap(allparcelinds)==Power(allparcelinds)) / numel(allparcelinds);

disp(['Overlap by area: ' num2str(nanmean(areaoverlap))])

templatefile = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc35175_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';
cifti_write_wHDR(output,templatefile,'Consensus_overlap_with_Power')