geo_distlim = 60;
euc_distlim = 40;
outname = '120_avg_corr_LR_smallwall_distregress_bspline_truesmooth';
smparm=.3; %(lower = stiffer)
nsmooth=10;

maskL = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/L.atlasroi_group_noproj.func.gii');maskL = maskL.cdata;%maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
maskR = gifti('/data/hcp-zfs/home/laumannt/120_parcellation/R.atlasroi_group_noproj.func.gii');maskR = maskR.cdata;%maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
% surfaceareaL = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR_surfaceareas.func.gii']);
% surfaceareaR = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR_surfaceareas.func.gii']);
% 
% surfacearea_cifti = zeros(66697,1);
% surfacearea_cifti(1:(nnz(maskL)+nnz(maskR))) = [surfaceareaL.cdata(logical(maskL)) ; surfaceareaR.cdata(logical(maskR))];

distmat = '/data/cn4/evan/Temp/distmat_smallwall_homo_surf_geodesic_vol_euc.mat';%'/data/hcp-zfs/home/laumannt/120_parcellation/modified_cifti_network/normalwall_distmat/distmat_homo_surf_geodesic_vol_euc.mat';%'/data/hcp-zfs/home/laumannt/120_parcellation/modified_cifti_network/normalwall_distmat/distmat_surf_geodesic_vol_euc.mat';%
corrfile = '/data/hcp-zfs/home/laumannt/120_parcellation/gradients_cifti_smallwall_wateredge/120/avg_corr_LR.func.gii';
outtemplate = '/data/hcp-zfs/home/laumannt/120_parcellation/gradients_cifti_smallwall_wateredge/120/avg_corr_LR.func.gii';

%%
load(distmat);
nrows = size(distmat,1);
withinhem_inds = logical(triu(ones(size(distmat_use)),1));
withinhem_inds(:,(nnz(maskL)+nnz(maskR)+1:end)) = 0;
withinhem_inds(1:nnz(maskL),(nnz(maskL)+1):(nnz(maskL)+nnz(maskR))) = 0;

distmat_col = distmat_use(logical(withinhem_inds));
clear distmat_use

corrmat = gifti(corrfile);
corrmat_col = corrmat.cdata(logical(withinhem_inds));
clear corrmat

% surfacearea_big = repmat(surfacearea_cifti,1,66697);
% surfacearea_col = surfacearea_big(logical(withinhem_inds));
% clear surfacearea_big
% [B, Bint, corrmat_col] = regress(corrmat_col,[surfacearea_col,ones(size(surfacearea_col))]);
% 
% binsize = [];
% indsub = [];
% fitted = zeros(size(corrmat_col));
% for i = 2:2:geo_distlim
%     bininds = find((distmat_col >= (i-2)) .* (distmat_col<i));
%     if ~isempty(bininds)
%         if isempty(binsize)
%             binsize = length(bininds);
%         end
%         randorder = randperm(length(bininds));
%         this_indsub = bininds(randorder(1:binsize));
%         indsub = [indsub ; this_indsub];
%         pp=bsmooth_evan(double(corrmat_col_noSA(this_indsub)),distmat_col(this_indsub),nsmooth,smparm);
%         thisfit=fnval(pp,distmat_col(bininds));
%         fitted(bininds) = thisfit;
%         thisresid = corrmat_col_noSA(bininds) - thisfit;
%         
%         figure;subplot(1,2,1);plot(distmat_col(this_indsub),corrmat_col_noSA(this_indsub),'.')
%         hold on
%         %subfitted = thisfit(this_indsub);
%         [sortdist, sorti] = sort(distmat_col(bininds),'ascend');
%         subplot(1,2,1);plot(sortdist,thisfit(sorti),'k-','LineWidth',3)
%         
%         subplot(1,2,2);plot(distmat_col(bininds),thisresid,'.')
%         
%     end
% end

%[B, Bint, corrmat_col_noSA] = regress(corrmat_col(indsub),[surfacearea_col(indsub),ones(size(indsub))]);

ind = find(distmat_col <= geo_distlim);
indsub = ind(1:1000:end);
%cftool(distmat_col(indsub),corrmat_col(indsub))
%curve = fit(distmat_col(indsub),corrmat_col(indsub),'power2');


pp=bsmooth_evan(double(corrmat_col(indsub)),double(distmat_col(indsub)),nsmooth,smparm);
fitted=fnval(pp,distmat_col);
fitted(distmat_col > geo_distlim) = 0;
withinhem_resid = corrmat_col - fitted;



figure;subplot(1,2,1);plot(distmat_col(indsub),corrmat_col(indsub),'.')
hold on
subfitted = fitted(indsub);
[sortdist, sorti] = sort(distmat_col(indsub),'ascend');
subplot(1,2,1);plot(sortdist,subfitted(sorti),'k-','LineWidth',3)

subplot(1,2,2);plot(distmat_col(indsub),withinhem_resid(indsub),'.')

% %%
% out = zeros(nrows);
% out(withinhem_inds) = withinhem_resid;
% 
% 
% out = out + out';
% 
% 
% cifti_write_wHDR(out,'/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.func.gii',outname,'dconn')


%%
load(distmat)
crosshem_inds = logical(zeros(size(distmat_use)));
crosshem_inds(1:nnz(maskL),(nnz(maskL)+1):(nnz(maskL)+nnz(maskR))) = 1;

[crosshem_replace(:,1), crosshem_replace(:,2)] = find((distmat_use==0) .* crosshem_inds);
crosshem_inds = logical(crosshem_inds .* distmat_use);

distmat_col = distmat_use(logical(crosshem_inds));
clear distmat_use
corrmat = gifti(corrfile);
corrmat_col = corrmat.cdata(logical(crosshem_inds));
clear corrmat
ind = find(distmat_col <= geo_distlim);
indsub = ind(1:1000:end);

% surfacearea_big = repmat(surfacearea_cifti,1,66697);
% surfacearea_col = surfacearea_big(logical(crosshem_inds));
% clear surfacearea_big
% [B, Bint, corrmat_col] = regress(corrmat_col,[surfacearea_col,ones(size(surfacearea_col))]);


% binsize = [];
% indsub = [];
% for i = 2:2:geo_distlim
%     bininds = find((distmat_col >= (i-2)) .* (distmat_col<i));
%     if ~isempty(bininds)
%         if isempty(binsize)
%             binsize = length(bininds);
%         end
%         randorder = randperm(length(bininds));
%         indsub = [indsub ; bininds(randorder(1:binsize))];
%     end
% end


pp=bsmooth_evan(double(corrmat_col(indsub)),double(distmat_col(indsub)),nsmooth,smparm);
fitted=fnval(pp,distmat_col);
fitted(distmat_col > geo_distlim) = 0;
crosshem_resid = corrmat_col - fitted;

figure;subplot(1,2,1);plot(distmat_col(indsub),corrmat_col(indsub),'.')
hold on
subfitted = fitted(indsub);
[sortdist, sorti] = sort(distmat_col(indsub),'ascend');
subplot(1,2,1);plot(sortdist,subfitted(sorti),'k-','LineWidth',3)

subplot(1,2,2);plot(distmat_col(indsub),crosshem_resid(indsub),'.')

%%
load(distmat)
volsurf_inds = logical(zeros(size(distmat_use)));
volsurf_inds(1:(nnz(maskL)+nnz(maskR)),(nnz(maskL)+nnz(maskR)+1):end) = 1;
distmat_col = distmat_use(logical(volsurf_inds));
clear distmat_use
corrmat = gifti(corrfile);
corrmat_col = corrmat.cdata(logical(volsurf_inds));
clear corrmat
ind = find(distmat_col <= euc_distlim);
indsub = ind(1:1000:end);



% binsize = [];
% indsub = [];
% for i = 2:2:euc_distlim
%     bininds = find((distmat_col >= (i-2)) .* (distmat_col<i));
%     if ~isempty(bininds)
%         if isempty(binsize)
%             binsize = length(bininds);
%         end
%         randorder = randperm(length(bininds));
%         indsub = [indsub ; bininds(randorder(1:binsize))];
%     end
% end

pp=bsmooth_evan(double(corrmat_col(indsub)),double(distmat_col(indsub)),nsmooth,smparm);
fitted=fnval(pp,distmat_col);
fitted(distmat_col > euc_distlim) = 0;
volsurf_resid = corrmat_col - fitted;

figure;subplot(1,2,1);plot(distmat_col(indsub),corrmat_col(indsub),'.')
hold on
subfitted = fitted(indsub);
[sortdist, sorti] = sort(distmat_col(indsub),'ascend');
subplot(1,2,1);plot(sortdist,subfitted(sorti),'k-','LineWidth',3)

subplot(1,2,2);plot(distmat_col(indsub),volsurf_resid(indsub),'.')

%%
load(distmat)
volvol_inds = logical(triu(ones(size(distmat_use)),1));
volvol_inds(1:(nnz(maskL)+nnz(maskR)),:) = 0;
distmat_col = distmat_use(logical(volvol_inds));
clear distmat_use
corrmat = gifti(corrfile);
corrmat_col = corrmat.cdata(logical(volvol_inds));
clear corrmat

ind = find(distmat_col <= euc_distlim);
indsub = ind(1:1000:end);



% binsize = [];
% indsub = [];
% for i = 2:2:euc_distlim
%     bininds = find((distmat_col >= (i-2)) .* (distmat_col<i));
%     if ~isempty(bininds)
%         if isempty(binsize)
%             binsize = length(bininds);
%         end
%         randorder = randperm(length(bininds));
%         indsub = [indsub ; bininds(randorder(1:binsize))];
%     end
% end

pp=bsmooth_evan(double(corrmat_col(indsub)),double(distmat_col(indsub)),nsmooth,smparm);
fitted=fnval(pp,distmat_col);
fitted(distmat_col > euc_distlim) = 0;
volvol_resid = corrmat_col - fitted;

figure;subplot(1,2,1);plot(distmat_col(indsub),corrmat_col(indsub),'.')
hold on
subfitted = fitted(indsub);
[sortdist, sorti] = sort(distmat_col(indsub),'ascend');
subplot(1,2,1);plot(sortdist,subfitted(sorti),'k-','LineWidth',3)

subplot(1,2,2);plot(distmat_col(indsub),volvol_resid(indsub),'.')

%%

out = zeros(size(crosshem_inds));
out(withinhem_inds) = withinhem_resid;
out(crosshem_inds) = crosshem_resid;
out(volsurf_inds) = volsurf_resid;
out(volvol_inds) = volvol_resid;

clear withinhem_inds crosshem_inds volsurf_inds volvol_inds withinhem_resid crosshem_resid volsurf_resid volvol_resid

load /data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat
for i = 1:size(crosshem_replace,1)
    nodeneighs = neighbors(crosshem_replace(i,2),2:7); nodeneighs(isnan(nodeneighs)) = [];
    out(crosshem_replace(i,1),crosshem_replace(i,2)) = mean(out(crosshem_replace(i,1),nodeneighs));
end

out = out + out';

for i = 1:(nnz(maskL)+nnz(maskR))
    nodeneighs = neighbors(i,2:7); nodeneighs(isnan(nodeneighs)) = [];
    out(i,i) = mean(out(i,nodeneighs));
end

for i = (nnz(maskL)+nnz(maskR)+1):size(out,1)
    out(i,i) = FisherTransform(1);
end

cifti_write_wHDR(out,outtemplate,outname,'dconn')


