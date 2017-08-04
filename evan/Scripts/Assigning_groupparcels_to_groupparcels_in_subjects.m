
hem = 'L';

subjectnums = [];

max_centroiddistance = 20;

datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';
iscifti = 2;

outputfolder = '/data/cn4/evan/RestingState/Ind_variability/Subjects/'; %location match data will be written to

parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_watershedmerge.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;
origparcels = parcels;

% baddata = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); baddata = baddata.cdata;
% baddata = baddata<750;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];

% parcelIDs_incrap = unique(parcels .* baddata);
% for badparcelID = parcelIDs_incrap'
%     parcelIDs(parcelIDs == badparcelID) = [];
% end

%load(['/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_allparcels/parcelcorrelpatterns_' hem '.mat'])
%parcelcorrelpatterns = squeeze(mean(parcelcorrelpatterns,1));
load(['/data/cn4/evan/RestingState/FC_Mapping_120/parcelcorrelpatterns_120_' hem '_watershedmerge.mat'])

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

[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

if ~isempty(subjectnums)
    subjects = subjects(subjectnums); subdata = subdata(subjectnums); tmasks = tmasks(subjectnums);
end


groupsimilarityoutput = zeros(length(mask),length(subjects));

%subject loop
for s = 1:length(subjects)
    
    
    disp(['Subject ' num2str(s)])
    
    
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{s} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{s});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    groupvsubcorrel = zeros(length(parcelIDs),1);
    
    for parcelnum = 1:length(parcelIDs)
        disp(['  Parcel number ' num2str(parcelnum)])
        
        parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + strcmp('R',hem)*ncortexLverts;
        
        parceltimecourse = mean(subtimecourse(parcelindices{parcelnum},:),1);
        subparcelcorrelpattern = paircorr_mod(parceltimecourse',subtimecourse');
        subparcelcorrelpattern(isnan(subparcelcorrelpattern)) = 0;
        
        indicestocorrelate = logical(ones(length(subparcelcorrelpattern),1)); indicestocorrelate(parcelindices{parcelnum}) = 0;
        
        similarity = paircorr_mod(subparcelcorrelpattern(indicestocorrelate)',parcelcorrelpatterns(parcelnum,indicestocorrelate)');
        
        groupsimilarityoutput(origparcels==parcelIDs(parcelnum),s) = similarity;
        
    end
    
end

save(gifti(single(groupsimilarityoutput)),[outputfolder 'Group_parcels_assignedtogroupparcels_insubs_similarity_' hem '.func.gii'])
    
