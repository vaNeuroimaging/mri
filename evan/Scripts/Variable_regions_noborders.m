hems = {'L','R'};

mindist = 10;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


for hemnum = 1:length(hems)
    hem = hems{hemnum};

    consensus = gifti(['120_LR_minsize400_recolored_manualconsensus_' hem '_cleaned.func.gii']); consensus = consensus.cdata;
    
    clustersize = gifti(['Variability_' hem '_consensus_dice_clustersize.func.gii']); clustersize = clustersize.cdata;
    
    IDs = gifti(['Variability_' hem '_consensus_dice.func.gii']); IDs = IDs.cdata;
    
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])
    
    inds = find(clustersize(:,2)>20);
    
    variable_noborders = zeros(32492,1);
    
    for ind = inds(:)'
        thisID = IDs(ind,2);
        if min(geo_distances(ind,(consensus==thisID))) > mindist
            variable_noborders(ind) = thisID;
        end
    end
    
    alternateIDs = unique(variable_noborders); alternateIDs(alternateIDs==0) = [];
    
    out = zeros(32492,0);
    
    for i = 1:length(alternateIDs)
        try
            data = metric_cluster(variable_noborders,alternateIDs(i)-.5,alternateIDs(i)+.5,30);
            out(:,(end+1):(end+size(data,2))) = data;
        catch
        end
    end
    
    save(gifti(single(out)),['Variable_regions_20_noborders_' hem '_clusters.func.gii'])
    
end