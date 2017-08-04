
hems = {'L','R'};

hemnum = 1;
    hem = hems{hemnum};
    
    
datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
iscifti = 2;
%hem = 'R';

%kdenval = .2;

outputfolder = ['/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_allparcels/'];
%parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_crossthresh_watershedmerge.func.gii'];
parcelfilename = ['/data/cn4/evan/ROIs/264_ROIs_' hem '.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;


system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample ' parcelfilename ' /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC ' parcelfilename(1:end-9) '_164.func.gii -largest'])

parcels_upsampled = gifti([parcelfilename(1:end-9) '_164.func.gii']); parcels_upsampled = parcels_upsampled.cdata;

baddata = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); baddata = baddata.cdata;
baddata = baddata<750;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0) =[];

totaltimeseries = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_LR_timeseries.func.gii';
totaltmask = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_total_tmask.txt';

% parcelIDs_incrap = unique(parcels .* baddata);
% for badparcelID = parcelIDs_incrap'
%     parcelIDs(parcelIDs == badparcelID) = [];
% end

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
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


for parcelnum = 1:length(parcelIDs)
    parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
    
end



[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{1} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{1});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    nsubcortverts = size(subtimecourse,1) - (ncortexLverts+ncortexRverts);

totaltimeseries = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_LR_timeseries.func.gii';
    totaltmask = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_total_tmask.txt'; 
consensusfile = '/data/cn4/evan/RestingState/Consensus/ConsensusMapvFinal.dtseries.nii';
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' consensusfile ' Temp.func.gii']);
    consensusdata = gifti('Temp.func.gii'); consensusdata = consensusdata.cdata;
    if iscifti==2
        tempconsensus = zeros(ncortexLverts+ncortexRverts+nsubcortverts,1);
        
        giftispace_temp = zeros(32492,1);
        ciftispace1_maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'); ciftispace1_maskL = ciftispace1_maskL.cdata;
        giftispace_temp(ciftispace1_maskL==0) = consensusdata(1:nnz(ciftispace1_maskL==0));
        tempconsensus(1:ncortexLverts) = giftispace_temp(logical(maskL.cdata));
        
        giftispace_temp = zeros(32492,1);
        ciftispace1_maskR = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'); ciftispace1_maskR = ciftispace1_maskR.cdata;
        giftispace_temp(ciftispace1_maskR==0) = consensusdata(nnz(ciftispace1_maskL==0)+1:nnz(ciftispace1_maskL==0)+nnz(ciftispace1_maskR==0));
        tempconsensus(ncortexLverts+1:ncortexLverts+ncortexRverts) = giftispace_temp(logical(maskR.cdata));
        
        tempconsensus(ncortexLverts+ncortexRverts+1:end) = consensusdata(nnz(ciftispace1_maskL==0)+nnz(ciftispace1_maskR==0)+1:end);
        
        consensusdata = tempconsensus;
        clear tempconsensus
    end
    
    
    networkIDs = unique(consensusdata);
    networkIDs(networkIDs==18) = [];
    networkIDs(networkIDs<=0) = [];
    
%     totaltimeseriesdata = gifti(totaltimeseries); totaltimeseriesdata = totaltimeseriesdata.cdata;
%     totaltmaskdata = load(totaltmask);
%     totaltimeseriesdata = totaltimeseriesdata(:,logical(totaltmaskdata));
%     totaltimeseriesdata(isnan(totaltimeseriesdata)) = 0;
    
    correlation = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_LR_timeseries_corr.func.gii'); correlation = correlation.cdata;
    correlation(isnan(correlation)) = 0;
    load /data/cn4/evan/RestingState/cifti_vertex_distance_matrix_thresh20mm.mat
    
    
    networkcorrelpattern = zeros(size(correlation,1),length(networkIDs));
    for networknum = 1:length(networkIDs)
        disp(num2str(networknum))
        
        networkverts = consensusdata==networkIDs(networknum);
        
        networkcorrelpattern(:,networknum) = sum((correlation(:,networkverts) .* matrix(:,networkverts)) , 2) ./ sum(matrix(:,networkverts) , 2);
        
    end
    %networkcorrelpattern(isnan(networkcorrelpattern)) = 0;
    
    save('/data/cn4/evan/RestingState/Consensus/Network_avgconnectivity_xdist20.mat','networkcorrelpattern','-v7.3')