function metric = find_center_water_parcel(watershedfile, hem)

outname = [watershedfile(1:end-9) '_centers.func.gii'];

water = gifti(watershedfile);
    water = water.cdata;

thisdir = pwd;

%HEMS = {'L';'R'};

%for hem = 1:2
    clear indpos
    surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
    %dir = '/data/cn4/laumannt/watershed_network/FCPROCESS_NEW/enhanced_water';
    %cd(dir)
    
    [phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    
    
    
    
    %Mask watersheds
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = mask.cdata;
    water(logical(mask)) = 0;
    
    
    waternum = unique(water);
    waternum(waternum==0) = [];
    
    metric = zeros(32492,1);    
    
    for w = 1:length(waternum)
        
        ind = find(water==waternum(w));
        
        %     meanphi = mean(phi(ind));
        %     meantheta = mean(theta(ind));
        %
        %     coord = [meanphi meantheta];
        %     sphere_coords = [phi(ind) theta(ind)];
        %
        %     rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        %
        %     dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        
        
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        
        
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
        
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        %indpos(w) = ind(indval);
        metric(ind(indval)) = waternum(w);
        
    end
    
    
    %metric(indpos) = 1;
    
    %Mask metric
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = mask.cdata;
    metric(logical(mask)) = 0;   
    indpos = find(metric);
    
    %Save out midthickness coordinates of centers
    midthick = gifti([surfdir '/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
    center_coords{hem} = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    save(gifti(single(metric)),outname)
    
    cd(thisdir)
    %dlmwrite([dir '/avgcorrofcorr_smooth2.55_allgrad_' HEMS{hem} '_smooth2.55_edge_avg_minima3_watershedmerged_parcel_center.txt'],indpos,'delimiter','\n')
    %save(['avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center.txt'],'indpos')
%end

%center_coords_both = [center_coords{1};center_coords{2}];
%quickroifile(center_coords_both,[dir '/avgcorrofcorr_smooth2.55_allgrad_LR_smooth2.55_edge_avg_minima3_watershedmerged_parcel_center.roi'])



    