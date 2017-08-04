bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

extend = 350;

outmetric = zeros(32492,1);

for vert = 1:12;
    
    outmetric(vert) = 1;
    
    for direction = 1:5;
        
        thisvert = vert;
        
        for i = 1:extend
        
            thisvert = neighbors(thisvert,direction+1);
            
            outmetric(thisvert) = 1;
            
        end
    end
end

save(gifti(single(outmetric)),'/data/cn4/evan/Temp/extend12.func.gii')
            