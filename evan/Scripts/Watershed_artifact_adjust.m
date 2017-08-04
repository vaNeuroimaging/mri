
factor = 1.05;

edgemapfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients_wateredge/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_wateredge_avg.func.gii';
edgemap = gifti(edgemapfile); edgemap = edgemap.cdata;

artifactlinefile = '/data/cn4/evan/Temp/Lines_between_first12.func.gii';
artifactlines = gifti(artifactlinefile); artifactlines = artifactlines.cdata;
artifactvertices = find(artifactlines);

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

outedgemap = edgemap;

lineratio = [];
borderlineratio = [];

for vert = artifactvertices';
    
    vertneighbors = neighbors(vert,2:7);
    vertneighbors(isnan(vertneighbors)) = [];
    
    nonlineneighbors = setdiff(vertneighbors,artifactvertices);
    
    newvalue = (edgemap(vert) + mean(edgemap(nonlineneighbors))) /2;
    
    lineratio(end+1) = edgemap(vert) / mean(edgemap(nonlineneighbors));
    
    outedgemap(vert) = newvalue;
    
    %outedgemap(vert) = edgemap(vert) *factor;
    
    for bordervert = nonlineneighbors
        
        bordervertneighbors = neighbors(bordervert,2:7);
        bordervertneighbors(isnan(bordervertneighbors)) = [];
        
        lineneighbors = intersect(bordervertneighbors,artifactvertices);
        
        %disp(length(lineneighbors))
        
        newvalue = (edgemap(bordervert) + mean(edgemap(lineneighbors))) / 2;
        %newvalue = mean([edgemap(bordervert) ; (edgemap(lineneighbors))]);
        
        borderlineratio(end+1) = edgemap(bordervert) / mean(edgemap(lineneighbors));
        
        outedgemap(bordervert) = newvalue;
        %outedgemap(bordervert) = edgemap(bordervert) /factor;
        
    end
    
end

save(gifti(single(outedgemap)),'/data/cn4/evan/Temp/AdjustedWatershedEdges.func.gii')
        
    