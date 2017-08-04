function spatial_corr_metric_Evan(metric1name,metric2name,neighdist,outputdir,outputname)

bufsize=16384;

a = dlmread(['/data/cn4/evan/fsaverage_LR32k/node_neighbors' num2str(neighdist) '.txt']);
maxneighbors = size(a,2);

outputstring = '[';
for i = 1:maxneighbors
    outputstring = [outputstring 'neighbors(:,' num2str(i) ') '];
end
outputstring = [outputstring(1:end-1) ']'];

eval([outputstring ' = textread([''/data/cn4/evan/fsaverage_LR32k/node_neighbors'' num2str(neighdist) ''.txt''],repmat(''%u '',[1 maxneighbors]),''delimiter'','' '',''bufsize'',bufsize,''emptyvalue'',NaN);']);

neighbors = neighbors+1;

metric1 = gifti(metric1name);
metric1 = metric1.cdata;

metric2 = gifti(metric2name);
metric2 = metric2.cdata;

corr_metric = zeros(size(metric1));


for i = 1:length(metric1)
    
    nodeneigh = neighbors(i,:);
    nodeneigh(isnan(nodeneigh)) = [];
    
    metric1_val = metric1(nodeneigh);
    metric2_val = metric2(nodeneigh);
    corr_metric(i) = paircorr_mod(metric1_val,metric2_val);
    
end

  save(gifti(single(corr_metric)),[outputdir '/' outputname '.func.gii']);
 
    
            

