vertexwise_regularized_fileL = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_vertexwise_infomap/geodist/120_Consensus_L_allthresh_recolored.func.gii';
parcelwise_regularized_fileL = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_infomap/watershed_Tk0002to005in00002_S1to1_surfxd30_INFMAP/120_L_minsize2_consensus_allthresh_recolor.func.gii';
parcels_fileL = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii';
parcelthreshs = [.002:.0002:.05];
parcelminthreshnum = 63;
vertthreshs = [.01:.002:.03];

vertexwise_regularized_fileR = '/data/cn4/evan/RestingState/FC_Mapping_120/120_R_vertexwise_infomap/eucdist/120_Consensus_R_allthresh.func.gii';
parcelwise_regularized_fileR = '/data/cn4/evan/RestingState/FC_Mapping_120/120_R_infomap/watershed_Tk0002to005in00002_S1to1_surfxd30_INFMAP/120_R_minsize2_consensus_allthresh_2.func.gii';
parcels_fileR = '/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.45.func.gii';
% parcelthreshs = [.002:.0002:.05];
% parcelminthreshnum = 53;
% vertthreshs = [.01:.002:.03];

%parcels = gifti(parcels_file); parcels = parcels.cdata;
%parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];

ciftitemplatefile = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');
maskR = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii');
template = gifti(ciftitemplatefile); template = template.cdata;


vertexwiseL = gifti(vertexwise_regularized_fileL); vertexwiseL = vertexwiseL.cdata;
parcelwiseL = gifti(parcelwise_regularized_fileL); parcelwiseL = parcelwiseL.cdata;

vertexwiseR = gifti(vertexwise_regularized_fileR); vertexwiseR = vertexwiseR.cdata;
parcelwiseR = gifti(parcelwise_regularized_fileR); parcelwiseR = parcelwiseR.cdata;

output = zeros(size(template,1),size(vertexwiseL,2));
output2 = zeros(size(template,1),size(vertexwiseL,2));
ciftibraininds = [1:(nnz(maskL.cdata==0) + nnz(maskR.cdata==0))];

realparcelthreshs = parcelthreshs(parcelminthreshnum:end);

for v = 1:size(vertexwiseL,2)
    for p = 1:size(parcelwiseL,2)
        parcelindsL = logical((parcelwiseL(:,p) > 0));% .* (vertexwise(:,v)>0));
        numoverlapL = nnz(parcelwiseL(parcelindsL,p)==vertexwiseL(parcelindsL,v));
        
        parcelindsR = logical((parcelwiseR(:,p) > 0));% .* (vertexwise(:,v)>0));
        numoverlapR = nnz(parcelwiseR(parcelindsR,p)==vertexwiseR(parcelindsR,v));
        
        pctoverlap(p) =  (numoverlapL + numoverlapR) / (nnz(parcelindsL) + nnz(parcelindsR));
    end
    
%     for p = 1:size(parcelwise,2)
%         parcelmatchpcts = zeros(length(parcelIDs),1);
%         for parcelnum = 1:length(parcelIDs);
%             parcelinds = logical(parcels==parcelIDs(parcelnum));
%             if mean(parcelwise(parcelinds,p)) > 0;
%                 parcelmatchpcts(parcelnum) = nnz(parcelwise(parcelinds,p)==vertexwise(parcelinds,v)) / nnz(parcelinds);
%             end
%         end
%         parcelmatchpcts(parcelmatchpcts==0) = [];
%             
%         pctoverlap(p) = mean(parcelmatchpcts); 
%     end
    
    
    [maxpct maxi] = max(pctoverlap);
    output(ciftibraininds,v) = [parcelwiseL(logical(maskL.cdata==0),maxi) ; parcelwiseR(logical(maskR.cdata==0),maxi)];
    disp(['Vertexwise thresh ' num2str(v) ' best overlap: parcel thresh ' num2str(maxi) ', overlap = ' num2str(maxpct)])
    %parcelwise(:,maxi) = 1000;
    matchingthreshnum = find(single(realparcelthreshs)==single(vertthreshs(v)));
    disp(['Vertexwise thresh ' num2str(v) ' same thresh: overlap = ' num2str(pctoverlap(matchingthreshnum))])
    if ~isempty(matchingthreshnum)
        output2(ciftibraininds,v) = [parcelwiseL(logical(maskL.cdata==0),matchingthreshnum) ; parcelwiseR(logical(maskR.cdata==0),matchingthreshnum)];
    end
end

cifti_write_wHDR(output,ciftitemplatefile,'120_Consensus_LR_allthresh_matchingparcelthreshs')
cifti_write_wHDR(output2,ciftitemplatefile,'120_Consensus_LR_allthresh_sameparcelthreshs')


    

