function varargout = bsmooth_evan(dep,indep,nsmooth,varargin)
if length(varargin)>0
    smparm=varargin{1};
else
    smparm=2e-4;
end
[~,I]=sort(indep);
%Sorted distances
dvs=indep(I);
%Sorted correlation values
Rvs=dep(I);
%Smooth with median and average filters
Ss=smooth(medfilt1(Rvs,nsmooth),nsmooth);
%Truncate ends
Sst=Ss(nsmooth:end-nsmooth);
dvst=dvs(nsmooth:end-nsmooth);
% Avoid repeated entries in input to csaps
[~,ia,~] = unique(dvst);
pp=csaps(dvst(ia),Sst(ia),smparm);
% Correct extrapolation function
pp=fnxtr(pp);
% % Evaluate spline system
% Rss=fnval(pp,dvs);
% % Project
% Rfin(I)=Rvs-Rss;
% 
% varargout{1}=Rfin;
% varargout{2}=dvs;
% varargout{3}=Rvs;
% varargout{4}=Rss;
varargout{1}=pp;
