geo_distlim = 60;
euc_distlim = 30;
outname = '120_avg_corr_LR_distregress3';

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
Linds = [1:32492]'; Linds = Linds(logical(maskL));
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
Rinds = [1:32492]'; Rinds = Rinds(logical(maskR));
%%
load /data/hcp-zfs/home/laumannt/120_parcellation/modified_cifti_network/normalwall_distmat/distmat_homo_surf_geodesic_vol_euc.mat
crosshem_inds = logical(zeros(size(distmat_use)));
crosshem_inds(1:29696,29697:59412) = 1;

[crosshem_replace(:,1), crosshem_replace(:,2)] = find((distmat_use==0) .* crosshem_inds);
crosshem_inds = logical(crosshem_inds .* distmat_use);

distmat_col = distmat_use(logical(crosshem_inds));
clear distmat_use
corrmat = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii');
corrmat_col = corrmat.cdata(logical(crosshem_inds));
clear corrmat
ind = find(distmat_col <= geo_distlim);
indsub = ind(1:1000:end);
%cftool(distmat_col(indsub),corrmat_col(indsub))
curve = fit(distmat_col(indsub),corrmat_col(indsub),'power2');

fitted = curve.a .* (distmat_col .^ curve.b) + curve.c;
fitted(distmat_col > geo_distlim) = 0;
crosshem_resid = corrmat_col - fitted;
figure;plot(distmat_col(1:10000:end),crosshem_resid(1:10000:end),'.')

%%
load /data/hcp-zfs/home/laumannt/120_parcellation/modified_cifti_network/normalwall_distmat/distmat_homo_surf_geodesic_vol_euc.mat
withinhem_inds = logical(triu(ones(size(distmat_use)),1));
withinhem_inds(:,59413:end) = 0;
withinhem_inds(1:29696,29697:59412) = 0;
distmat_col = distmat_use(logical(withinhem_inds));
clear distmat_use
corrmat = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii');
corrmat_col = corrmat.cdata(logical(withinhem_inds));
clear corrmat
ind = find(distmat_col <= geo_distlim);
indsub = ind(1:1000:end);
%cftool(distmat_col(indsub),corrmat_col(indsub))
curve = fit(distmat_col(indsub),corrmat_col(indsub),'power2');

fitted = curve.a .* (distmat_col .^ curve.b) + curve.c;
fitted(distmat_col > geo_distlim) = 0;
withinhem_resid = corrmat_col - fitted;
figure;plot(distmat_col(1:10000:end),withinhem_resid(1:10000:end),'.')



%%
load /data/hcp-zfs/home/laumannt/120_parcellation/modified_cifti_network/normalwall_distmat/distmat_homo_surf_geodesic_vol_euc.mat
volsurf_inds = logical(zeros(size(distmat_use)));
volsurf_inds(1:59412,59413:end) = 1;
distmat_col = distmat_use(logical(volsurf_inds));
clear distmat_use
corrmat = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii');
corrmat_col = corrmat.cdata(logical(volsurf_inds));
clear corrmat
ind = find(distmat_col <= euc_distlim);
indsub = ind(1:100:end);
%cftool(distmat_col(indsub),corrmat_col(indsub))
curve = fit(distmat_col(indsub),corrmat_col(indsub),'power2');

fitted = curve.a .* (distmat_col .^ curve.b) + curve.c;
fitted(distmat_col > euc_distlim) = 0;
volsurf_resid = corrmat_col - fitted;
figure;plot(distmat_col(1:1000:end),volsurf_resid(1:1000:end),'.')

%%
load /data/hcp-zfs/home/laumannt/120_parcellation/modified_cifti_network/normalwall_distmat/distmat_homo_surf_geodesic_vol_euc.mat
volvol_inds = logical(triu(ones(size(distmat_use)),1));
volvol_inds(1:59412,:) = 0;
distmat_col = distmat_use(logical(volvol_inds));
clear distmat_use
corrmat = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii');
corrmat_col = corrmat.cdata(logical(volvol_inds));
clear corrmat
ind = find(distmat_col <= euc_distlim);
indsub = ind(1:100:end);
%cftool(distmat_col(indsub),corrmat_col(indsub))
curve = fit(distmat_col(indsub),corrmat_col(indsub),'power2');

fitted = curve.a .* (distmat_col .^ curve.b) + curve.c;
fitted(distmat_col > euc_distlim) = 0;
volvol_resid = corrmat_col - fitted;
figure;plot(distmat_col(1:1000:end),volvol_resid(1:1000:end),'.')

%%

out = zeros(66697);
out(withinhem_inds) = withinhem_resid;
out(crosshem_inds) = crosshem_resid;
out(volsurf_inds) = volsurf_resid;
out(volvol_inds) = volvol_resid;


load /data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat
for i = 1:size(crosshem_replace,1)
    nodeneighs = neighbors(crosshem_replace(i,2),2:7); nodeneighs(isnan(nodeneighs)) = [];
    out(crosshem_replace(i,1),crosshem_replace(i,2)) = mean(out(crosshem_replace(i,1),nodeneighs));
end

out = out + out';

for i = 1:59412
    nodeneighs = neighbors(i,2:7); nodeneighs(isnan(nodeneighs)) = [];
    out(i,i) = mean(out(i,nodeneighs));
end

for i = 59413:size(out,1)
    out(i,i) = FisherTransform(1);
end

cifti_write_wHDR(out,'/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii',outname,'dconn')


