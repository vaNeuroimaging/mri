rois = textread('/data/cn4/evan/ROIs/264_cortical_coords.txt');

hems = {'L','R'};
hemsign = [-1 1];



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
    
    verts{hemnum} = [];
    
    surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    midthick = gifti([surfdir '/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
    
    medialmaskdata = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']);
    medialmaskdata = medialmaskdata.cdata;
    
    for roinum = 1:size(rois,1)
        coords = rois(roinum,:);
        
        if sign(coords(1))==hemsign(hemnum);
            
            rep_coord = repmat(coords, [size(midthick.vertices,1) 1]);
            
            dist_coord = sum((midthick.vertices-rep_coord).^2,2).^(1/2);
            [ign index] = min(dist_coord);
            
            verts{hemnum} = [verts{hemnum}; index];
        end
    end
    
    dlmwrite(['/data/cn4/evan/ROIs/264_surfverts_' hem '.txt'],verts{hemnum})
    system(['wb_command -surface-geodesic-rois ' surfdir '/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii 10 /data/cn4/evan/ROIs/264_surfverts_' hem '.txt /data/cn4/evan/ROIs/264_surfverts_' hem '.func.gii -overlap-logic CLOSEST'])
    
    label = gifti(['/data/cn4/evan/ROIs/264_surfverts_' hem '.func.gii']); label = label.cdata;
    
    output = zeros(32492,1);
    for i = 1:size(label,2);
        output(logical(label(:,i))) = i;
    end
    save(gifti(single(output)),['/data/cn4/evan/ROIs/264_surfvert_ROIs_' hem '.func.gii'])
    
%     label = zeros(32492,1);
%     label(verts{hemnum}) = verts{hemnum};
%     
%     stillexpanding = 1;
%     while stillexpanding==1
%         stillexpanding = 0;
%         borderverts = find((label==0) .* (medialmaskdata==0));
%         borderverts = borderverts(randperm(length(borderverts)));
%         for vert = borderverts'
%             vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs)) = [];
%             neighlabels = label(vertneighs); neighlabels(neighlabels==0) = [];
%             if nnz(unique(neighlabels)) == 1
%                 stillexpanding = 1;
%                 label(vert) = neighlabels(1);
%             end
%         end
%     end
%             
%        
%     save(gifti(single(label)),['/data/cn4/evan/ROIs/264_surfvert_parcels_' hem '.func.gii'])
    
end





