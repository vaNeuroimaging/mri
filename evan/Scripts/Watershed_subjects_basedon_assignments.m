tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';
[subjects tmasks] = textread(tmasklist,'%s%s');
s=1;

thresholds = [.3 : .01 : .6];

evalc(['!wb_command -cifti-convert -to-gifti-ext /data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subjects{s} '/' subjects{s} '_timeseries_normalwall.dtseries.nii Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    delete('Temp.func.gii*');
        tmask = load(tmasks{s});
        subtimecourse = subtimecourse(:,logical(tmask));
    subtimecourse(isnan(subtimecourse)) = 0;

hem = 'L';

max_centroiddistance = 20;

min_correl = .05;

expandifneeded = 1;



parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;
origparcels = parcels;

randomizedparcels = gifti(['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45_randomized.func.gii']);
randomizedparcels = randomizedparcels.cdata;

parcelcorrelpatternfile = ['/data/cn4/evan/RestingState/FC_Mapping_120/Group_parcelcorrelpatterns_' hem '.mat'];

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];

load(parcelcorrelpatternfile)
clear neighbors
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);


load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);
surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
midthick = gifti([surfdir '/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
group_parcel_centroid_vert = zeros(length(parcelIDs),1);

for parcelnum = 1:length(parcelIDs)
    ind = find(origparcels==parcelIDs(parcelnum));
    randomizedIDs(parcelnum) = mean(randomizedparcels(ind));
    parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
    
    
    meanX = mean(sphere.vertices(ind,1));
    meanY = mean(sphere.vertices(ind,2));
    meanZ = mean(sphere.vertices(ind,3));
    
    coord = [meanX meanY meanZ];
    sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
    
    rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
    
    dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
    [y indval] = min(dist_coord);
    group_parcel_centroid_vert(parcelnum) = ind(indval);
end

groupsimilarities = zeros(length(thresholds),length(parcelIDs));

groupassignoutput = zeros(length(mask),length(thresholds));
groupsimilarityoutput = zeros(length(mask),length(thresholds));
subassignoutput = zeros(length(mask),length(thresholds));
subsimilarityoutput = zeros(length(mask),length(thresholds));


for threshnum = 1:length(thresholds)
    threshperc = thresholds(threshnum);
    
    
    groupassignmap = zeros(nnz(mask),1);
groupsimilaritymap = zeros(nnz(mask),1);
groupdistancemap = zeros(nnz(mask),1);
subassignmap = zeros(nnz(mask),1);
subsimilaritymap = zeros(nnz(mask),1);
subdistancemap = zeros(nnz(mask),1);
    
    watershed_algorithm_merge_andthresh_perc(['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subjects{s} '/avgcorrofcorr_allgrad_L_smooth2.55_wateredge_thresh0.01_avg.func.gii'],'./','Temp',hem,threshperc);
    subparcels = gifti(['Tempwatershedmerge_' num2str(threshperc) '.func.gii']); subparcels = subparcels.cdata;
    origsubparcels = subparcels;
    subparcels = subparcels(logical(mask));
    subparcelIDs = unique(subparcels); subparcelIDs(subparcelIDs==0)=[];
    
    groupvsubcorrel = zeros(length(parcelIDs),length(subparcelIDs));
    
    sub_parcel_centroid_vert = zeros(length(subparcelIDs),1);
    
    groupassignmentIDs = zeros(length(parcelIDs),1);
    subassignmentIDs = zeros(length(subparcelIDs),1);
    
    for subparcelnum = 1:length(subparcelIDs)
        string{subparcelnum} = ['  Parcel number ' num2str(subparcelnum) ' out of ' num2str(length(subparcelIDs))];
        if subparcelnum==1; fprintf('%s',string{subparcelnum}); else fprintf([repmat('\b',1,length(string{subparcelnum-1})) '%s'],string{subparcelnum}); end
        
        
        subparcelindices{subparcelnum} = find(subparcels==subparcelIDs(subparcelnum)) + (strcmp(hem,'R') * ncortexLverts);
        
        parceltimecourse = mean(subtimecourse(subparcelindices{subparcelnum},:),1);
        subparcelcorrelpattern = paircorr_mod(parceltimecourse',subtimecourse');
        subparcelcorrelpattern(isnan(subparcelcorrelpattern)) = 0;
        
        ind = find(origsubparcels==subparcelIDs(subparcelnum));
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
        
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        sub_parcel_centroid_vert(subparcelnum) = ind(indval);
        
        thisparcel_distances = geo_distances(group_parcel_centroid_vert,sub_parcel_centroid_vert(subparcelnum));
        
        
        for parcelnum = 1:length(parcelIDs)
            indicestoavoid = union(find(geo_distances(:,group_parcel_centroid_vert(parcelnum))<20), find(geo_distances(:,sub_parcel_centroid_vert(subparcelnum))<20));
            indicestocorrelate = logical(ones(length(subparcelcorrelpattern),1)); indicestocorrelate(indicestoavoid) = 0;
            groupvsubcorrel(parcelnum,subparcelnum) = paircorr_mod(parcelcorrelpatterns(parcelnum,indicestocorrelate)',subparcelcorrelpattern(indicestocorrelate)');
        end
        
        
        thisparcel_groupvsubcorrel = groupvsubcorrel(:,subparcelnum);
        thisparcel_groupvsubcorrel(thisparcel_distances>max_centroiddistance) = -Inf;
        thisparcel_groupvsubcorrel(thisparcel_groupvsubcorrel<min_correl) = -Inf;
        
        [correl groupparcel] = max(thisparcel_groupvsubcorrel);
        if expandifneeded
            expanded_centroiddistance = max_centroiddistance;
            while correl==-Inf
                expanded_centroiddistance = expanded_centroiddistance + 5;
                thisparcel_groupvsubcorrel = groupvsubcorrel(:,subparcelnum);
                thisparcel_groupvsubcorrel(thisparcel_distances>expanded_centroiddistance) = -Inf;
                %if expanded_centroiddistance < 50
                    thisparcel_groupvsubcorrel(thisparcel_groupvsubcorrel<min_correl) = -Inf;
                %end
                [correl groupparcel] = max(thisparcel_groupvsubcorrel);
            end
        elseif correl==-Inf
            error(['No parcels within ' num2str(max_centroiddistance) ' vertices!'])
        end
        subassignment(subparcelnum) = groupparcel;

        if groupassignmentIDs(groupparcel)
            subassignmentIDs(subparcelnum) = randomizedIDs(groupparcel);
            groupassignment{groupparcel} = [groupassignment{groupparcel} subparcelnum];
        else
            subassignmentIDs(subparcelnum) = randomizedIDs(groupparcel);
            groupassignmentIDs(groupparcel) = randomizedIDs(groupparcel);
            groupassignment{groupparcel} = subparcelnum;
        end
        
    end
    
    for parcelnum = 1:length(parcelIDs)
    
    if (groupassignmentIDs(parcelnum) > 0) %&& (onevtwo_pctoverlap(parcelnum,groupassignment(parcelnum)) > 0)
        
        groupassignmap(parcels==parcelIDs(parcelnum)) = groupassignmentIDs(parcelnum);
        
            groupsimilaritymap(parcels==parcelIDs(parcelnum)) = mean(groupvsubcorrel(parcelnum,groupassignment{groupparcel}));
            groupsimilarities(threshnum,parcelnum) = mean(groupvsubcorrel(parcelnum,groupassignment{groupparcel}));
        
        
    else
        
        groupassignmap(parcels==parcelIDs(parcelnum)) = -1;
        groupsimilaritymap(parcels==parcelIDs(parcelnum)) = -1;
        groupdistancemap(parcels==parcelIDs(parcelnum)) = -1;
    end
end

for parcelnum = 1:length(subparcelIDs)
    
        
        subassignmap(subparcels==subparcelIDs(parcelnum)) = subassignmentIDs(parcelnum);
        
        subsimilaritymap(subparcels==subparcelIDs(parcelnum)) = mean(groupvsubcorrel(subassignment(parcelnum),parcelnum));
    
        

end

   subassignmap_temp = zeros(32492,1);
   subassignmap_temp(logical(mask)) = subassignmap;
   subassignmap = subassignmap_temp;
   groupassignmap_temp = zeros(32492,1);
   groupassignmap_temp(logical(mask)) = groupassignmap;
   groupassignmap = groupassignmap_temp;
   groupassignoutput(:,threshnum) = groupassignmap;
   groupsimilarityoutput(logical(mask),threshnum) = groupsimilaritymap;
   subassignoutput(:,threshnum) = subassignmap;
   subsimilarityoutput(logical(mask),threshnum) = subsimilaritymap;
   
   overall_assignment_similarity(threshnum) = mean(groupsimilarities(threshnum,:));
end

[ign maxi] = max(overall_assignment_similarity);

disp(['Best threshold: ' num2str(thresholds(maxi))])

copyfile(['Temp_watershedmerge_' num2str(thresholds(maxi)) '.func.gii'],[subjects{s} '_bestassignment_watersheds.func.gii'])
delete('Temp*.func.gii')

save(gifti(single(subassignoutput(:,maxi))),[subjects{s} '_subtogroup_assignment.func.gii']);
    