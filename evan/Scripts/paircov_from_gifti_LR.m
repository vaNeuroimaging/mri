function [cov_corr] = paircov_from_gifti_LR(giftiname,hem)
%PAIRCORR Computes pairwise Pearson's linear correlation coefficient with
% optional significance. Returns r, a p1-by-p2 matrix containing the
% pairwise correlation coefficient between each pair of columns in the
% n-by-p1 and n-by-p2 matrices a and b. r is calculated as the dot
% product between two vectors divided by the product of their magnitudes.
% If a second output argument is provided, like so:
% [r p] = paircorr(a,b)
% then p is the two-tailed significance.
% Added single input functionality TOL, 04/01/12.

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = ~(mask.cdata);

masknameL = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']; maskL = gifti(masknameL); ncortLverts = nnz(~maskL.cdata);

a = gifti(giftiname); a = a.cdata(:,(1:nnz(mask))+(strcmp(hem,'R') * ncortLverts));

cov_corr = cov(a);
