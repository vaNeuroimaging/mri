datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
iscifti = 2;
hem = 'L';

variableparcelIDs = [15024];

assignmentsfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_latOccip/watershed_bothhem_Tk02to02in001_S1to1_BI_INFMAP/rawassn_minsize10.txt';
assignmentcolumn = 1;
communityIDs = [1];


outputfolder = ['/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_latOccip/'];

parcels = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/120_L_crossthresh_watershedmerge.func.gii'); parcels = parcels.cdata;

baddata = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); baddata = baddata.cdata;
baddata = baddata<800;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];

% parcelIDs_incrap = unique(parcels .* baddata);
% for badparcelID = parcelIDs_incrap'
%     parcelIDs(parcelIDs == badparcelID) = [];
% end

% subject_correlations_ofparcels = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/Subject_avgcorrelation_to_120_L_crossthresh_watershedmerge.func.gii');
% subject_correlations_ofparcels = subject_correlations_ofparcels.cdata;
 for parcelnum = 1:length(parcelIDs)
%     variableparcelindex(parcelnum) = mean(subject_correlations_ofparcels(parcels==parcelIDs(parcelnum))) < .55;
    variableparcelindex(parcelnum) = any(variableparcelIDs==parcelIDs(parcelnum));
 end


%origparcels = gifti(parcelfilename); origparcels = origparcels.cdata;



%%

if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');maskL.cdata = ~logical(maskL.cdata); ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii');maskR.cdata = ~logical(maskR.cdata);ncortexRverts = nnz(maskR.cdata);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    parcels = origparcels;
    ncortexLverts = length(parcels);
    ncortexRverts = length(parcels);
end


for parcelnum = 1:length(variableparcelIDs)
    parcelindices{parcelnum} = find(parcels==variableparcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
end

assignments = load(assignmentsfile); assignments = assignments(:,assignmentcolumn);

[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

for communityID = communityIDs;
    communitysubjectinds = find(assignments==communityID);
    communitysubjects = subjects(communitysubjectinds);
    
    for communitysubnum = 1:length(communitysubjects)
    
    subnum = communitysubjectinds(communitysubnum);
    disp(['Subject ' communitysubjects{communitysubnum}])
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{subnum});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    parceltimecourses = zeros(length(variableparcelIDs),size(subtimecourse,2));
    
    for parcelnum = 1:length(variableparcelIDs)
        parceltimecourses(parcelnum,:) = mean(subtimecourse(parcelindices{parcelnum},:),1);
    end
    
    meanparcelconnectivity(communitysubnum,:) = paircorr_mod(mean(parceltimecourses,1)',subtimecourse');
        
    end

    community_parcelconnectivity = mean(meanparcelconnectivity,1);
    
    outputL = zeros(length(mask),1);
    outputL(logical(maskL.cdata)) = community_parcelconnectivity(1:ncortexLverts);
    save(gifti(single(outputL)),[outputfolder '/Community' num2str(communityID) '_connectivityL.func.gii']);
    
    
    outputR = zeros(length(mask),1);
    outputR(logical(maskR.cdata)) = community_parcelconnectivity(ncortexLverts+1:ncortexLverts+ncortexRverts);
    save(gifti(single(outputR)),[outputfolder '/Community' num2str(communityID) '_connectivityR.func.gii']);
    
    clear meanparcelconnectivity

end