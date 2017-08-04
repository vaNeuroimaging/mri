max_centroiddistance = 20;

min_correlval = 0;%.25;

hem = 'L';

iscifti = 2;

outputfolder = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/32_C1/'; %location watersheds and match data will be written to


parcelfilename1 = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/32_C1/C1_L_2.0_watershedmerge.func.gii';

timecoursefilename1 = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_lininterp_222/gradients_cifti_voxelclean_erode_concat_wateredge/C1_all64/allsubs_LR_timeseries.dtseries.nii';

tmaskfilename1 = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_lininterp_222/gradients_cifti_voxelclean_erode_concat_wateredge/C1_all64/allsubs_total_tmask.txt';


parcelfilename2 = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/32_C2/C2_L_2.0_watershedmerge.func.gii';

timecoursefilename2 = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_lininterp_222/gradients_cifti_voxelclean_erode_concat_wateredge/C2_all64/allsubs_LR_timeseries.dtseries.nii';

tmaskfilename2 = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_lininterp_222/gradients_cifti_voxelclean_erode_concat_wateredge/C2_all64/allsubs_total_tmask.txt';


parcels1 = gifti(parcelfilename1); parcels1 = parcels1.cdata;
origparcels1 = parcels1;
parcelIDs1 = unique(parcels1);
parcelIDs1(parcelIDs1==0)=[];

parcels2 = gifti(parcelfilename2); parcels2 = parcels2.cdata;
origparcels2 = parcels2;
parcelIDs2 = unique(parcels2); parcelIDs2(parcelIDs2==0)=[];

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels1 = parcels1(logical(mask));
    parcels2 = parcels2(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels1 = parcels1(logical(mask));
    parcels2 = parcels2(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    parcels1 = parcels_upsampled;
    ncortexLverts = length(parcels1);
end

%load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);

surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);

disp('Calculating parcel connectivities for map 1')

evalc(['!wb_command -cifti-convert -to-gifti-ext ' timecoursefilename1 ' Temp1.func.gii']);
timecourse1 = gifti('Temp1.func.gii'); timecourse1 = timecourse1.cdata;

    tmask = load(tmaskfilename1);
    timecourse1 = timecourse1(:,logical(tmask));
timecourse1(isnan(timecourse1)) = 0;


parcel1_centroid_vert = zeros(length(parcelIDs1),1);
for parcelnum = 1:length(parcelIDs1)
    %disp(['  Parcel number ' num2str(parcelnum)])
    
    parcelindices1{parcelnum} = find(parcels1==parcelIDs1(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
     
     parceltimecourse1(parcelnum,:) = mean(timecourse1(parcelindices1{parcelnum},:),1);
    
%     ind = find(origparcels1==parcelIDs1(parcelnum));
%     meanX = mean(sphere.vertices(ind,1));
%     meanY = mean(sphere.vertices(ind,2));
%     meanZ = mean(sphere.vertices(ind,3));
%     
%     coord = [meanX meanY meanZ];
%     sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
%     
%     rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
%     
%     dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
%     [y indval] = min(dist_coord);
%     parcel1_centroid_vert(parcelnum) = ind(indval);
end

parcelcorrelpattern1 = paircorr_mod(parceltimecourse1',timecourse1');
parcelcorrelpattern1(isnan(parcelcorrelpattern1)) = 0;


disp('Calculating parcel connectivities for map 2')

evalc(['!wb_command -cifti-convert -to-gifti-ext ' timecoursefilename2 ' Temp2.func.gii']);
timecourse2 = gifti('Temp2.func.gii'); timecourse2 = timecourse2.cdata;

    tmask = load(tmaskfilename2);
    timecourse2 = timecourse2(:,logical(tmask));
timecourse2(isnan(timecourse2)) = 0;


parcel2_centroid_vert = zeros(length(parcelIDs2),1);
for parcelnum = 1:length(parcelIDs2)
    %disp(['  Parcel number ' num2str(parcelnum)])
    
    parcelindices2{parcelnum} = find(parcels2==parcelIDs2(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
     
     parceltimecourse2(parcelnum,:) = mean(timecourse2(parcelindices2{parcelnum},:),1);
    
%     ind = find(origparcels2==parcelIDs2(parcelnum));
%     meanX = mean(sphere.vertices(ind,1));
%     meanY = mean(sphere.vertices(ind,2));
%     meanZ = mean(sphere.vertices(ind,3));
%     
%     coord = [meanX meanY meanZ];
%     sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
%     
%     rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
%     
%     dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
%     [y indval] = min(dist_coord);
%     parcel2_centroid_vert(parcelnum) = ind(indval);
end

onevtwo_pctoverlap = zeros(length(parcelIDs1),length(parcelIDs2));
for i = 1:length(parcelIDs1)
    for j = 1:length(parcelIDs2)
        numoverlapverts = length(intersect(parcelindices1{i},parcelindices2{j}));
        %onevtwo_pctoverlap(i,j) = mean([numoverlapverts/length(parcelindices1{i}),numoverlapverts/length(parcelindices2{j})]);
        onevtwo_pctoverlap(i,j) = numoverlapverts;
    end
end


 parcelcorrelpattern2 = paircorr_mod(parceltimecourse2',timecourse2');
 parcelcorrelpattern2(isnan(parcelcorrelpattern2)) = 0;
 
 disp('Conducting assignments and saving outputs')
 
 onevtwocorrel = paircorr_mod(parcelcorrelpattern1',parcelcorrelpattern2');

%onevtwo_distances = geo_distances(parcel1_centroid_vert,parcel2_centroid_vert);

onevtwocost = ((max(max(onevtwo_pctoverlap))-onevtwo_pctoverlap) / (max(max(onevtwo_pctoverlap)))) + (1-(onevtwocorrel+1)/2);
%onevtwocost(onevtwo_pctoverlap==0) = Inf;
%onevtwocost = (1-(onevtwocorrel+1)/2);
%onevtwocost(onevtwo_distances > max_centroiddistance) = Inf;
%onevtwocost(onevtwocorrel < min_correlval) = Inf;

%calculate the assignments
[assignment, cost] = munkres(onevtwocost);

randomizedIDs = randperm(length(parcelIDs1));

oneassignmap = zeros(nnz(mask),1);
onesimilaritymap = zeros(nnz(mask),1);
%onedistancemap = zeros(nnz(mask),1);
twoassignmap = zeros(nnz(mask),1);
twosimilaritymap = zeros(nnz(mask),1);
%twodistancemap = zeros(nnz(mask),1);

for parcelnum = 1:length(parcelIDs1)
    
    if (assignment(parcelnum) > 0) && (onevtwo_pctoverlap(parcelnum,assignment(parcelnum)) > 0)
        
        oneassignmap(parcels1==parcelIDs1(parcelnum)) = randomizedIDs(parcelnum);
        onesimilaritymap(parcels1==parcelIDs1(parcelnum)) = onevtwocorrel(parcelnum,assignment(parcelnum));
%        onedistancemap(parcels1==parcelIDs1(parcelnum)) = onevtwo_distances(parcelnum,assignment(parcelnum));
        
        twoassignmap(parcels2==parcelIDs2(assignment(parcelnum))) = randomizedIDs(parcelnum);
        twosimilaritymap(parcels2==parcelIDs2(assignment(parcelnum))) = onevtwocorrel(parcelnum,assignment(parcelnum));
%        twodistancemap(parcels2==parcelIDs2(assignment(parcelnum))) = onevtwo_distances(parcelnum,assignment(parcelnum));
        
    else
        
        oneassignmap(parcels1==parcelIDs1(parcelnum)) = -1;
        onesimilaritymap(parcels1==parcelIDs1(parcelnum)) = -1;
%        onedistancemap(parcels1==parcelIDs1(parcelnum)) = -1;
        
    end
end

for parcelnum = 1:length(parcelIDs2)
    
    if ~any(twoassignmap(parcels2==parcelIDs2(parcelnum)))
        
        twoassignmap(parcels2==parcelIDs2(parcelnum)) = -1;
        twosimilaritymap(parcels2==parcelIDs2(parcelnum)) = -1;
%        twodistancemap(parcels2==parcelIDs2(parcelnum)) = -1;
        
    end
end



assgnoutput1 = zeros(length(mask),1); assgnoutput1(logical(mask),:) = oneassignmap; origassgnoutput1 = assgnoutput1;
assgnoutput2 = zeros(length(mask),1); assgnoutput2(logical(mask),:) = twoassignmap; origassgnoutput2 = assgnoutput2;

simoutput1 = zeros(length(mask),1); simoutput1(logical(mask),:) = onesimilaritymap;
simoutput2 = zeros(length(mask),1); simoutput2(logical(mask),:) = twosimilaritymap;


for vertnum = 1:size(assgnoutput1,1)
    vertneighs = neighbors(vertnum,2:7);
    vertneighs(isnan(vertneighs)) = [];
    if origassgnoutput1(vertnum)==-1 && ~any(origassgnoutput1(vertneighs)==0)
        assgnoutput1(vertnum)=0;
        simoutput1(vertnum) = 0;
    end
    if origassgnoutput2(vertnum)==-1 && ~any(origassgnoutput2(vertneighs)==0)
        assgnoutput2(vertnum)=0;
        simoutput2(vertnum) = 0;
    end
end


save(gifti(single(assgnoutput1)),[outputfolder 'Map1_parcels_assignedtosubs_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])
save(gifti(single(assgnoutput2)),[outputfolder 'Map2_parcels_assignedtogroup_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])
save(gifti(single(simoutput1)),[outputfolder 'Map1_parcels_assignedtosubs_similarity_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])
save(gifti(single(simoutput2)),[outputfolder 'Map2_parcels_assignedtogroup_similarity_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])

%output = zeros(length(mask),1); output(logical(mask),:) = onesimilaritymap;


%output = zeros(length(mask),1); output(logical(mask),:) = onedistancemap;
%save(gifti(single(output)),[outputfolder 'Map1_parcels_assignedtosubs_distance_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])


%output = zeros(length(mask),1); output(logical(mask),:) = twosimilaritymap;


%output = zeros(length(mask),1); output(logical(mask),:) = twodistancemap;
%save(gifti(single(output)),[outputfolder 'Map2_parcels_assignedtogroup_distance_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])


