iscifti = 2;

outputfolder = ['/data/cn4/evan/RestingState/FC_Mapping_120/'];
hems = {'L','R'};
for hemnum = 1:2
    hem = hems{hemnum};
    parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_crossthresh_watershedmerge.func.gii'];
    %parcelfilename = ['/data/cn4/evan/ROIs/264_ROIs_' hem '.func.gii'];
    parcels = gifti(parcelfilename); parcels = parcels.cdata;
    origparcels = parcels;
    
    parcelIDs = unique(parcels);
    parcelIDs(parcelIDs==0) =[];
    
    totaltimeseries = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_LR_timeseries.func.gii';
    totaltmask = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_total_tmask.txt';
    
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
    
    timeseries = gifti(totaltimeseries); timeseries = timeseries.cdata;
    tmask = load(totaltmask);
    timeseries = timeseries(:,logical(tmask));
    
    load(['parcelcorrelpatterns_' hem '.mat'])
    
    meanoutput = zeros(length(mask),1);
    stdoutput = zeros(length(mask),1);
    
    for parcelnum = 1:length(parcelIDs)
        disp(['Parcel ' num2str(parcelnum)])
        parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
        
        parcelcorrelpattern = paircorr_mod(timeseries',mean(timeseries(parcelindices{parcelnum},:),1)');
        parcelcorrelpattern(isnan(parcelcorrelpattern)) = 0;
        
        subjectcorrelations = paircorr_mod(squeeze(parcelcorrelpatterns(:,parcelnum,:))',parcelcorrelpattern);
        
        meancorrelation = mean(subjectcorrelations);
        meanoutput(origparcels==parcelIDs(parcelnum)) = meancorrelation;
        stdcorrelation = std(subjectcorrelations);
        stdoutput(origparcels==parcelIDs(parcelnum)) = stdcorrelation;
        
    end
    
    save(gifti(single(meanoutput)),['Parcel_mean_correlation_wgroup_' hem '.func.gii'])
    save(gifti(single(stdoutput)),['Parcel_std_correlation_wgroup_' hem '.func.gii'])
    
end