function [bold,tempbetas] = demean_detrend(bold,tmask)
if ~exist('tmask')
    tmask = true(1,size(bold,2));
end
[vox,ts] = size(bold);
linreg=[repmat(1,[ts 1]) linspace(0,1,ts)'];
tempboldcell=num2cell(bold(:,logical(tmask))',1);
linregcell=repmat({linreg(logical(tmask),:)},[1 vox]);
tempbetas = cellfun(@mldivide,linregcell,tempboldcell,'uniformoutput',0);
tempbetas=cell2mat(tempbetas);
tempbetas=tempbetas';
tempintvals=tempbetas*linreg';
bold=bold-tempintvals;

end