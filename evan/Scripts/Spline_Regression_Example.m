%To tune for real data:
%Adjust number of subsamples and smoothing parameter
nsub=10000;
smparm=.3; %(lower = stiffer)

%this approach will struggle with a R vs d relationship that
%is very flat but sharply increases at very low distances
%spline regressing the near (<20mm?) and far separately may be optimal

%number of ROIs
n=1000; 
%a distance matrix
d=10*rand(n)+1;
d(~~eye(size(d)))=0;
%a (z) correlation matrix with some 1/R effect
R=.05*randn(n)+1./(d);
R(~~eye(size(R)))=inf;
% length of initial smooth/medfilt1/edge truncation to remove outliers
nsmooth=10;%fix(numel(R)./25);

%remove diagonal
d_nd=d(~eye(size(d)));
R_nd=R(~eye(size(R)));

%random subsampling of elements
subsample=unique(randi(numel(R_nd),nsub,1));

%compute the spline system (pp)
[Rout,ds,Rs,Rss,pp]=bsmooth(R_nd(subsample),d_nd(subsample),nsmooth,smparm);
%Outputs:
%ds:   sorted distance matrix
%Rs:   correlation matrix (sorted by distance)
%Rss:  regression line (sorted by distance)
% Rss is ds evaluated by the spline system pp
%Rout: detrended input (unsorted back to original order of input)

%evaluate spline system on full data structure
R_nd_spline=fnval(pp,d_nd);
R_nd_detrend=R_nd-R_nd_spline;

%reconstruct matrix
R_detrend=zeros(size(R));
R_detrend(~eye(size(R)))=R_nd_detrend;

figure(8675309),clf
subplot 211,hold on
plot(d(:),R(:),'.')
plot(d_nd(:),R_nd_spline(:),'r.')
xlim([min(nonzeros(d)) max(nonzeros(d))])

subplot 212
plot(d(:),R_detrend(:),'.')
xlim([min(nonzeros(d)) max(nonzeros(d))])