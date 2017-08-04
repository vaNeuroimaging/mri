%% Find watershed parcel centers

center_coords_both = [];

%all_water_corrmat is the correlation matrix (only used to get the size)
xdistancemat = ones(size(all_water_corrmat));

distanceexclusion = 30;

parceldir = '/data/cn5/selfRegulation/V4Process_nosmooth/gradients_120_108_combined_subsurf_nosmooth/';
parcelnames = {'120_108_combined_L_watershedmerge_0.35_tweaked.func.gii','120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'};%,'Poldrome_R_crossthresh_watershedmerge.func.gii'};


for hem = 1:length(HEMS)
    clear indpos
    surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    midthick = gifti([surfdir '/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    
    [phi theta r] = cart2sph(midthick.vertices(:,1), midthick.vertices(:,2),midthick.vertices(:,3));
    
    thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']); 
    thismedialwall = ~thismedialwall.cdata;
    
    parcel = gifti([parceldir '/' parcelnames{hem}]);
    parcel = parcel.cdata;    
    
    parcelnum = unique(parcel);
    parcelnum(parcelnum==0) = [];
    
    for w = 1:length(parcelnum)
        
        ind = find(parcel==parcelnum(w));        
        
        meanX = mean(midthick.vertices(ind,1));
        meanY = mean(midthick.vertices(ind,2));
        meanZ = mean(midthick.vertices(ind,3));
        
        
        coord = [meanX meanY meanZ];
        midthick_coords = [midthick.vertices(ind,1) midthick.vertices(ind,2) midthick.vertices(ind,3)];
        
        rep_coord = repmat(coord, [size(midthick_coords,1) 1]);
        
        dist_coord = sum((midthick_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        indpos(w) = ind(indval);
        
    end
    
    metric = zeros(32492,1);    
    metric(indpos) = 1;
    
    %Save out midthickness coordinates of centers
    midthick = gifti([surfdir '/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    
    if distanceexclusion > 0
        load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' HEMS{hem} '.mat'])
        hemxdistancemat = geo_distances(indpos,indpos);
        clear geo_distances indpos
        hemxdistancemat = hemxdistancemat > distanceexclusion;
        
        if hem==1
            heminds = 1:size(hemxdistancemat,1);
        else
            heminds = (size(xdistancemat,1) -size(hemxdistancemat,1) +1) : size(xdistancemat,1);
        end
        xdistancemat(heminds,heminds) = hemxdistancemat;
    end
    
    center_coords_both = [center_coords_both ; center_coords];
end