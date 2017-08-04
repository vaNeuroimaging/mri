

hems = {'L','R'};
for hemnum = 1:length(hems)
    origROIs = textread('/data/cn4/evan/ROIs/PowerROIs.txt');
    
    hem = hems{hemnum};
    
    xdirection = (hemnum*2)-3;
    
    ROIs = origROIs(((origROIs(:,1).*xdirection)>0),:);
    if hemnum==1
        ROIs = [ROIs ; origROIs((origROIs(:,1)==0),:)];
    end
    
    surf = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
    surf_verts = surf.vertices;
    
    alldistances = distance(ROIs',surf_verts');
    
    [ign minvertices] = min(alldistances,[],2);
    
    minvertices = minvertices-1;
    
    dlmwrite(['/data/cn4/evan/ROIs/264_verts_' hem '.txt'],minvertices)
    
    system(['wb_command -surface-geodesic-rois /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii 5 /data/cn4/evan/ROIs/264_verts_' hem '.txt /data/cn4/evan/ROIs/264_ROIs_' hem '.func.gii -overlap-logic CLOSEST'])
    
    surf_ROIs = gifti(['/data/cn4/evan/ROIs/264_ROIs_' hem '.func.gii']); surf_ROIs = surf_ROIs.cdata;
    surf_parcels = zeros(size(surf_ROIs,1),1);
    
    for i = 1:size(surf_ROIs,2)
        surf_parcels(logical(surf_ROIs(:,i))) = i;
    end
    
    save(gifti(single(surf_parcels)),['/data/cn4/evan/ROIs/264_ROIs_' hem '.func.gii']);
end




    