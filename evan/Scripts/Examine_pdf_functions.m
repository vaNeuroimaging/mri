origdata = gifti('Poldrome_R.func.gii'); origdata = origdata.cdata;
mask = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'); mask = ~mask.cdata;
data = double(origdata(logical(mask)));

bins = 100;
bincentervals = ([1/bins:1/bins:1]-(1/2/bins));
binedgevals = ([0:1/bins:1]);

% funcname = 'ex-gaussian';
% exgaussparams = simple_egfit(data);
%fittedcurve = exgausspdf(exgaussparams(1),exgaussparams(2),exgaussparams(3),bincentervals);


funcname = 'gamma';
params = fitdist(data,funcname);
if params.NumParams==3
    fittedcurve = pdf(funcname,bincentervals,params.Params(1),params.Params(2),params.Params(3));
elseif params.NumParams==2
    fittedcurve = pdf(funcname,bincentervals,params.Params(1),params.Params(2));
else
    fittedcurve = pdf(funcname,bincentervals,params.Params(1));
end
counts = histc(data,binedgevals);
counts(end-1) = counts(end-1)+counts(end); counts = counts(1:end-1);
bar(bincentervals,counts,'histc')
hold on
plot(bincentervals,counts,'g-')
plot(bincentervals,fittedcurve/100*nnz(mask),'r-')
hold off
disp([funcname ' function explains ' num2str((corr2(fittedcurve/100*nnz(mask),counts')^2)*100) '% of variance'])