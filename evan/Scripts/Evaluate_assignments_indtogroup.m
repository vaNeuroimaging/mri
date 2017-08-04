datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';

outputfolder = '/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/';

%distanceexclusion = 30;

hems = {'L','R'};

subjectnums = [2:3];

[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

if ~isempty(subjectnums)
    subjects = subjects(subjectnums); subdata = subdata(subjectnums); tmasks = tmasks(subjectnums);
end

assignedconnectivitysimilarity = [];
connectivitysimilarity = [];

counter = 1;

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
    %maskR = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii');ncortexRverts = nnz(maskR.cdata==0);

    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);
    geo_distances = geo_distances(:,logical(mask));
    
    sphere = gifti(['/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
    
    g2sassignments = gifti(['/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Group_parcels_assignedtosubs_' hem '.func.gii']);
    g2sassignments = g2sassignments.cdata(logical(mask),:);
    s2gassignments = gifti(['/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Subject_parcels_assignedtogroup_' hem '.func.gii']);
    s2gassignments = s2gassignments.cdata(logical(mask),:);
    
    parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45.func.gii'];
    parcels = gifti(parcelfilename); 
    origparcels = parcels.cdata;
    parcels = parcels.cdata(logical(mask));
    parcels(g2sassignments(:,1)==0) = 0;
    parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
    
    assigned_correlpatterns = zeros(66697,length(subjects),length(parcelIDs));
    correlpatterns = zeros(66697,length(subjects),length(parcelIDs));
    
    for s = 1:length(subjects)
        
        evalc(['!wb_command -cifti-convert -to-gifti-ext /data/cn4/evan/RestingState/Ind_variability/Subjects/' subjects{s} '/' subjects{s} '_timeseries_normalwall.dtseries.nii Temp.func.gii']);
        subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
        delete('Temp.func.gii*');
        if ~isempty(tmasklist)
            tmask = load(tmasks{s});
            subtimecourse = subtimecourse(:,logical(tmask));
        end
        subtimecourse(isnan(subtimecourse)) = 0;
        
        for parcelnum = 1:length(parcelIDs)
            
            string{counter} = ['Subject ' num2str(s) ', ' hem 'hem, parcel number ' num2str(parcelnum)];
            if counter==1; fprintf('%s',string{counter}); else fprintf([repmat('\b',1,length(string{counter-1})) '%s'],string{counter}); end
            counter = counter+1;
            
            assignID = mean(g2sassignments(parcels==parcelIDs(parcelnum),s));
            indices = find(s2gassignments(:,s)==assignID) + (ncortexLverts*strcmp('R',hem));
            assigned_parceltimecourse = mean(subtimecourse(indices,:),1);
            this_correlpattern = paircorr_mod(subtimecourse',assigned_parceltimecourse');
            if any(~isnan(this_correlpattern));
                this_correlpattern(isnan(this_correlpattern)) = 0;
            end
            assigned_correlpatterns(:,s,parcelnum) = this_correlpattern;
            
            indices = find(parcels==(parcelIDs(parcelnum))) + (ncortexLverts*strcmp('R',hem));
            parceltimecourse = mean(subtimecourse(indices,:),1);
            this_correlpattern = paircorr_mod(subtimecourse',parceltimecourse');
            if any(~isnan(this_correlpattern));
                this_correlpattern(isnan(this_correlpattern)) = 0;
            end
            correlpatterns(:,s,parcelnum) = this_correlpattern;
        end
    end
    
    
    mean_assigned_correlpatterns = zeros(32492,length(parcelIDs));
    mean_correlpatterns = zeros(32492,length(parcelIDs));
    
%     mean_assigned_correlpatterns = mean(assigned_correlpatterns((1:nnz(mask)) + (ncortexLverts*strcmp('R',hem)),:,:),2);
%     mean_assigned_correlpatterns = squeeze(mean_assigned_correlpatterns);
%     mean_correlpatterns = mean(correlpatterns((1:nnz(mask)) + (ncortexLverts*strcmp('R',hem)),:,:),2);
%     mean_correlpatterns = squeeze(mean_correlpatterns);
    
    thishem_assignedconnectivitysimilarity = zeros(length(parcelIDs),1);
    thishem_assignedconnectivitysimilarity_map = zeros(size(mask));
    
    thishem_connectivitysimilarity = zeros(length(parcelIDs),1);
    thishem_connectivitysimilarity_map = zeros(size(mask));
    
    for parcelnum = 1:length(parcelIDs)
        
        subswiththisparcel = ~isnan(sum(assigned_correlpatterns(:,:,parcelnum),1));
        mean_assigned_correlpatterns(logical(mask),parcelnum) = mean(assigned_correlpatterns((1:nnz(mask)) + (ncortexLverts*strcmp('R',hem)),subswiththisparcel,parcelnum),2);
        mean_correlpatterns(logical(mask),parcelnum) = mean(correlpatterns((1:nnz(mask)) + (ncortexLverts*strcmp('R',hem)),subswiththisparcel,parcelnum),2);
        
        
        ind = find(origparcels==parcelIDs(parcelnum));
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
    
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
    
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
    
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        group_parcel_centroid_vert = ind(indval);
        
        vertices_tocorrelate = (geo_distances(group_parcel_centroid_vert,:) > 20);
        
        crosscorrel = paircorr_mod(squeeze(assigned_correlpatterns(vertices_tocorrelate,:,parcelnum)));
        crosscorrel(isnan(crosscorrel)) = 0;
        offdiagvals = crosscorrel(~diag(logical(diag(crosscorrel))));
        avg_crosscorrel = mean(offdiagvals(~isnan(offdiagvals)));
        thishem_assignedconnectivitysimilarity(parcelnum) = avg_crosscorrel;
        thishem_assignedconnectivitysimilarity_map(origparcels==parcelIDs(parcelnum)) = avg_crosscorrel;
        
        crosscorrel = paircorr_mod(squeeze(correlpatterns(vertices_tocorrelate,:,parcelnum)));
        offdiagvals = crosscorrel(~diag(logical(diag(crosscorrel))));
        avg_crosscorrel = mean(offdiagvals(~isnan(offdiagvals)));
        thishem_connectivitysimilarity(parcelnum) = avg_crosscorrel;
        thishem_connectivitysimilarity_map(origparcels==parcelIDs(parcelnum)) = avg_crosscorrel;
%         if parcelnum==12
%             1;
%         end
    end
    
    assignedconnectivitysimilarity = [assignedconnectivitysimilarity; thishem_assignedconnectivitysimilarity];
    connectivitysimilarity = [connectivitysimilarity; thishem_connectivitysimilarity];
    
    save(gifti(single(thishem_assignedconnectivitysimilarity_map)),[outputfolder '/Assignedparcels_connectivity_similarity_' hem '.func.gii'])
    save(gifti(single(thishem_connectivitysimilarity_map)),[outputfolder '/Groupparcels_connectivity_similarity_' hem '.func.gii'])
    
    save(gifti(single(mean_assigned_correlpatterns)),[outputfolder '/Mean_assignedparcels_connectivity_' hem '.func.gii'])
    save(gifti(single(mean_correlpatterns)),[outputfolder '/Mean_groupparcels_connectivity_' hem '.func.gii'])
    
end
disp(' ')

[H,P,CI,STATS] = ttest2(assignedconnectivitysimilarity,connectivitysimilarity)
            
        
        
    
    