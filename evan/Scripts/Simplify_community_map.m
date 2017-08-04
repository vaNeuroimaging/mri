function Simplify_community_map(metricname,minsize,hem,outline)
%Simplify_community_map(metricname,minsize,hem,[outline])

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%load the gifti
% if isnumeric(metricname)
%     metric = metricname;
% else
metric = gifti(metricname);
metric = metric.cdata;
%end

clusters = communities_discrete_clusters(metric);

simplified = zeros(size(metric));

for cluster=1:size(clusters,2);
    if nnz(clusters(:,cluster)) > minsize;
        simplified(logical(clusters(:,cluster))) = clusters(logical(clusters(:,cluster)),cluster);
    end
end

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
medialinds = find(mask);

baddataname = ['/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/RotParcellation/Baddata_bigcluster_' hem '.func.gii'];
baddata = gifti(baddataname); baddata = baddata.cdata;
baddatainds = find(baddata);

badmedialinds = union(medialinds,baddatainds);

zeroinds = find(simplified==0);
zeroinds = setdiff(zeroinds,badmedialinds);
numzeros = nnz(zeroinds);
while ~isempty(zeroinds)
    tempsimplified = simplified;
    for vert = zeroinds(:)'
        vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs)) = [];
        neighvals = simplified(vertneighs); neighvals(neighvals==0) = [];
        if ~isempty(neighvals)
            tempsimplified(vert) = mode(neighvals);
        end
    end
    simplified = tempsimplified;
    zeroinds = find(simplified==0);
    zeroinds = setdiff(zeroinds,badmedialinds);
    if nnz(zeroinds)==numzeros;
        zeroinds = [];
    end
    numzeros = nnz(zeroinds);
end

save(gifti(single(simplified)),[metricname(1:end-9) '_simplified.func.gii'])

if exist('outline','var') && (outline>0)
    borders = zeros(size(simplified));
    
    for i = 1:length(borders)
        indneighs = neighbors(i,2:end); indneighs(isnan(indneighs)) = [];
        neighvals = simplified(indneighs); neighvals(neighvals<1) = [];
        if any(neighvals~=simplified(i))
            
            borders(i) = 1;
        end
    end
    save(gifti(single(borders)),[metricname(1:end-9) '_simplified_outline.func.gii'])
end