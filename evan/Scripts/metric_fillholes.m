function metricout = metric_fillholes(metric)

%function metric_fillholes(metricnamein,metricnameout)

% metric = gifti(metricnamein);
% metric = metric.cdata;
 metricout = metric;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

for map = 1:size(metric,2)
    for node = 1:size(metric,1)
        
        nodeneighs = neighbors(node,2:7);
        nodeneighs(isnan(nodeneighs)) = [];
        
        nodeval = metric(node,map);
        neighvals = metric(nodeneighs,map);
        
        if length(find(neighvals==nodeval)) < 2 && nodeval~=15;
           uniqueneighvals = unique(neighvals);
           maxnumber = 0;
           for i = 1:length(uniqueneighvals)
               if length(find(neighvals==uniqueneighvals(i))) > maxnumber
                   maxnumber = length(find(neighvals==uniqueneighvals(i)));
                   newnodeval = uniqueneighvals(i);
               end
           end
           metricout(node,map) = newnodeval;
        end
    end
end

%save(gifti(single(metricout)),metricnameout);
            