hem = 'L';

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);

outmetric = zeros(32492,1);

vertpairs = [1:12];

outmetric(vertpairs) = 1;

for verti = vertpairs;
    for vertj = vertpairs
        
        if verti > vertj;
            
            thisline = zeros(32492,1);
            
            disp(['Vertex ' num2str(verti) ' vs ' num2str(vertj)])
            
            vertj_spherecoords = sphere.vertices(vertj,:);
            
            thisvert = verti;
            
            distancetraveled = 0;
            
            while thisvert ~= vertj
                
                thisvertneighs = neighbors(thisvert,2:7);
                thisvertneighs(isnan(thisvertneighs)) = [];
                
                thisneighspherecoords = sphere.vertices(thisvertneighs',:);
                thisneigh_coorddiffs = thisneighspherecoords - repmat(vertj_spherecoords,[length(thisvertneighs) 1]);
                thisneigh_distaway = sqrt((thisneigh_coorddiffs(:,1).^2) + (thisneigh_coorddiffs(:,2).^2) + (thisneigh_coorddiffs(:,3).^2));
                [mindist minind] = min(thisneigh_distaway);
                thisvert = thisvertneighs(minind);
                
                distancetraveled = distancetraveled + 1;
                
                thisline(thisvert) = 1;
                
            end
            
            if distancetraveled < 70;
                
                outmetric = outmetric + thisline;
                
            end
            
        end
    end
end

save(gifti(single(outmetric)),'/data/cn4/evan/Temp/Lines_between_first12.func.gii')
                
                

