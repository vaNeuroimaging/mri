
hems = {'L','R'};
for hemnum = 1:2
hem = hems{hemnum};

subjectnums = [];

max_centroiddistance = 20;

min_correl = .05;
buffer = .05;

expandifneeded = 1;

datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';
iscifti = 1;

outputfolder = '/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Matching_closestgood/'; %location match data will be written to

parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;
origparcels = parcels;

randomizedparcels = gifti(['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45_randomized.func.gii']);
randomizedparcels = randomizedparcels.cdata;

parcelcorrelpatternfile = ['/data/cn4/evan/RestingState/FC_Mapping_120/Group_parcelcorrelpatterns_' hem '.mat'];

% baddata = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); baddata = baddata.cdata;
% baddata = baddata<750;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];

% parcelIDs_incrap = unique(parcels .* baddata);
% for badparcelID = parcelIDs_incrap'
%     parcelIDs(parcelIDs == badparcelID) = [];
% end

load(parcelcorrelpatternfile)
clear neighbors
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    parcels = parcels_upsampled;
    ncortexLverts = length(parcels);
end

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
%group_center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
%clear indpos

[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

if ~isempty(subjectnums)
    subjects = subjects(subjectnums); subdata = subdata(subjectnums); tmasks = tmasks(subjectnums);
end


%randomizedIDs = randperm(length(parcelIDs));

groupassignoutput = zeros(length(mask),length(subjects));
groupsimilarityoutput = zeros(length(mask),length(subjects));
groupdistanceoutput = zeros(length(mask),length(subjects));
subassignoutput = zeros(length(mask),length(subjects));
subsimilarityoutput = zeros(length(mask),length(subjects));
subdistanceoutput = zeros(length(mask),length(subjects));

%subject loop
for s = 1:length(subjects)
    
    groupassignmap = zeros(nnz(mask),1);
groupsimilaritymap = zeros(nnz(mask),1);
groupdistancemap = zeros(nnz(mask),1);
subassignmap = zeros(nnz(mask),1);
subsimilaritymap = zeros(nnz(mask),1);
subdistancemap = zeros(nnz(mask),1);



    
    disp(['Subject ' num2str(s)])
    
    subparcelfilename = ['/data/cn4/evan/RestingState/Ind_variability/Subjects/' subjects{s} '/' subjects{s} '_' hem '_watershedmerge_0.45.func.gii'];
    subparcels = gifti(subparcelfilename); subparcels = subparcels.cdata;
    origsubparcels = subparcels;
    subparcels = subparcels(logical(mask));
    subparcelIDs = unique(subparcels); subparcelIDs(subparcelIDs==0)=[];
    
    %randomizedIDs = randperm(length(parcelIDs)+length(subparcelIDs));

    groupassignmentIDs = zeros(length(parcelIDs),1);
    subassignmentIDs = zeros(length(subparcelIDs),1);
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext /data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subjects{s} '/' subjects{s} '_timeseries_normalwall.dtseries.nii Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    delete('Temp.func.gii*');
    if ~isempty(tmasklist)
        tmask = load(tmasks{s});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    groupvsubcorrel = zeros(length(parcelIDs),length(subparcelIDs));
    
    sub_parcel_centroid_vert = zeros(length(subparcelIDs),1);
    
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
        
        maxcorrel = max(thisparcel_groupvsubcorrel);
        goodmatches = find(thisparcel_groupvsubcorrel >= (maxcorrel-buffer));
        [ign closesti] = min(thisparcel_distances(goodmatches));
        groupparcel = goodmatches(closesti);
        
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
    disp(' ')
    
%     for parcelnum = 1:length(parcelIDs)
%         
%         if ~groupassignmentIDs(parcelnum)
%             thisgroupparcel_correls = groupvsubcorrel(parcelnum,:);
%             thisgroupparcel_correls(geo_distances(group_parcel_centroid_vert(parcelnum),sub_parcel_centroid_vert) > 20) = -Inf;
%             [ign bestmatchsubparcel] = max(thisgroupparcel_correls);
%             groupassignmentIDs(parcelnum) = subassignmentIDs(bestmatchsubparcel);
%         end
%     end
    
    





for parcelnum = 1:length(parcelIDs)
    
    if (groupassignmentIDs(parcelnum) > 0) %&& (onevtwo_pctoverlap(parcelnum,groupassignment(parcelnum)) > 0)
        
        groupassignmap(parcels==parcelIDs(parcelnum)) = groupassignmentIDs(parcelnum);
        
            groupsimilaritymap(parcels==parcelIDs(parcelnum)) = mean(groupvsubcorrel(parcelnum,groupassignment{groupparcel}));
            %groupdistancemap(parcels==parcelIDs(parcelnum)) = mean(group_v_sub_distance(parcelnum,(subassignment==parcelnum)));
        
        
    else
        
        groupassignmap(parcels==parcelIDs(parcelnum)) = -1;
        groupsimilaritymap(parcels==parcelIDs(parcelnum)) = -1;
        groupdistancemap(parcels==parcelIDs(parcelnum)) = -1;
    end
end

for parcelnum = 1:length(subparcelIDs)
    
    %if subassignment(parcelnum) > 0  %&& (onevtwo_pctoverlap(twoassignment(parcelnum),parcelnum) > 0)
        
        subassignmap(subparcels==subparcelIDs(parcelnum)) = subassignmentIDs(parcelnum);
        
        subsimilaritymap(subparcels==subparcelIDs(parcelnum)) = mean(groupvsubcorrel(subassignment(parcelnum),parcelnum));
        %subdistancemap(subparcels==subparcelIDs(parcelnum)) = mean(group_v_sub_distance((groupassignment==parcelnum),parcelnum));
        
        
    %else
%        
%         subassignmap(subparcels==subparcelIDs(parcelnum)) = -1;
%         subsimilaritymap(subparcels==subparcelIDs(parcelnum)) = -1;
%         subdistancemap(subparcels==subparcelIDs(parcelnum)) = -1;
%         
%     end
end

    
    
   subassignmap_temp = zeros(32492,1);
   subassignmap_temp(logical(mask)) = subassignmap;
   subassignmap = subassignmap_temp;
   groupassignmap_temp = zeros(32492,1);
   groupassignmap_temp(logical(mask)) = groupassignmap;
   groupassignmap = groupassignmap_temp;

% done = 0;
% iter = 0;
% while done==0
%     iter = iter+1;
%     %disp(['Iteration ' num2str(iter)])
%     done = 1;
% 
%     subassignments = unique(subassignmap); subassignments(subassignments==0) = [];
%     allvals = [1:max(subassignments)];
%     for i = 1:length(subassignments)
%         vertices = find(subassignmap==subassignments(i));
%         neighbors = find(sum((geo_distances(vertices,:)<10),1));
%         neighbors = setdiff(neighbors,vertices);
%         neighborvals = unique(subassignmap(neighbors)); neighborvals(neighborvals==0) = []; 
%         if any(abs(neighborvals-subassignments(i))<15)
%             inside_neighborvalrange = [];
%             for thisneighval = neighborvals(:)'
%                 inside_neighborvalrange = unique([inside_neighborvalrange [thisneighval-15 : thisneighval+30]]);
%             end
%             remainingvals = setdiff(allvals,[subassignments(:)' inside_neighborvalrange]);
%             newval = remainingvals(randi(length(remainingvals),1));
%             subassignmap(subassignmap==subassignments(i)) = newval;
%             groupassignmap(groupassignmap==subassignments(i)) = newval;
%             done = 0;
%         end
%     end
% end


    
    

     
    clear subparcelindices
    
    groupassignoutput(:,s) = groupassignmap;
    groupsimilarityoutput(logical(mask),s) = groupsimilaritymap;
    %groupdistanceoutput(logical(mask),s) = groupdistancemap;
    subassignoutput(:,s) = subassignmap;
    subsimilarityoutput(logical(mask),s) = subsimilaritymap;
    %subdistanceoutput(logical(mask),s) = subdistancemap;
    

    

    
end







save(gifti(single(groupassignoutput)),[outputfolder 'Group_parcels_matchtosubs_closestgood_' hem '.func.gii'])

save(gifti(single(groupsimilarityoutput)),[outputfolder 'Group_parcels_matchtosubs_closestgood_similarity_' hem '.func.gii'])

%save(gifti(single(groupdistanceoutput)),[outputfolder 'Group_parcels_assignedtosubs_distance_' hem '.func.gii'])

save(gifti(single(subassignoutput)),[outputfolder 'Subject_parcels_matchtogroup_closestgood_' hem '.func.gii'])

save(gifti(single(subsimilarityoutput)),[outputfolder 'Subject_parcels_matchtogroup_closestgood_similarity_' hem '.func.gii'])

%save(gifti(single(subdistanceoutput)),[outputfolder 'Subject_parcels_assignedtogroup_distance_' hem '.func.gii'])

parcelexistancemap = zeros(32492,1);
parceloverlapmap = zeros(32492,length(parcelIDs));

for parcelnum = 1:length(parcelIDs)
    
    for s = 1:length(subjects)
        matchID = mean(groupassignoutput(origparcels==parcelIDs(parcelnum),s));
        thissubparcelindices = find(subassignoutput(:,s)==matchID);
        parceloverlapmap(thissubparcelindices,parcelnum) = parceloverlapmap(thissubparcelindices,parcelnum) + 1;
        
        if ~isempty(thissubparcelindices)
            parcelexistancemap(origparcels==parcelIDs(parcelnum)) = parcelexistancemap(origparcels==parcelIDs(parcelnum))+1;
        end
    end
    parceloverlapmap(:,parcelnum) = parceloverlapmap(:,parcelnum) ./ mean(parcelexistancemap(origparcels==parcelIDs(parcelnum)));
end

save(gifti(single(parcelexistancemap)),[outputfolder 'Group_parcels_matchnumber_closestgood_' hem '.func.gii']);
save(gifti(single(parceloverlapmap)),[outputfolder 'Sub_parcels_overlap_closestgood_' hem '.func.gii']);
    
    
end