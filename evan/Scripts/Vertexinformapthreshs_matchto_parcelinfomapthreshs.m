% vertexwise_regularized_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_vertexwise_infomap/geodist/120_Consensus_L_allthresh_recolored.func.gii';
% parcelwise_regularized_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_infomap/watershed_Tk0002to005in00002_S1to1_surfxd30_INFMAP/120_L_minsize2_consensus_allthresh_recolor.func.gii';
% parcels_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii';

vertexwise_regularized_file = '/data/cn4/evan/Temp/120_Consensus_LR_allthresh4.dtseries.nii';
parcelwise_regularized_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_LR_infomap/watershed_Tk0002to005in00002_S1to1_surfxd30_INFMAP/120_LR_minsize4_consensus_allthresh2.dtseries.nii';
parcels_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_LR_wateredgethresh_watershedmerge_0.45_gooddata.dtseries.nii';

iscifti = 1;

parcelthreshs = [.002:.0002:.05];
parcelminthreshnum = 42;
vertthreshs = [.01:.002:.03];

maxthreshdist = .005;

% vertexwise_regularized_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_R_vertexwise_infomap/eucdist/120_Consensus_R_allthresh.func.gii';
% parcelwise_regularized_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_R_infomap/watershed_Tk0002to005in00002_S1to1_surfxd30_INFMAP/120_R_minsize2_consensus_allthresh.func.gii';
% parcels_file = '/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.45.func.gii';
% parcelthreshs = [.002:.0002:.05];
% parcelminthreshnum = 53;
% vertthreshs = [.01:.002:.03];


if iscifti
    parcels = cifti_read(parcels_file);

vertexwise = cifti_read(vertexwise_regularized_file);
parcelwise = cifti_read(parcelwise_regularized_file); 
else

parcels = gifti(parcels_file); parcels = parcels.cdata;

vertexwise = gifti(vertexwise_regularized_file); vertexwise = vertexwise.cdata;
parcelwise = gifti(parcelwise_regularized_file); parcelwise = parcelwise.cdata;
end

parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];

output = zeros(size(vertexwise));
output2 = zeros(size(vertexwise));
output3 = zeros(size(vertexwise));
overlap = zeros(size(parcelwise));
overlapoutput = zeros(size(vertexwise));

realparcelthreshs = parcelthreshs(parcelminthreshnum:end);
vert_pctoverlap = zeros(size(vertexwise,2),1);
for v = 1:size(vertexwise,2)
    for p = 1:size(parcelwise,2)
        parcelinds = logical((parcelwise(:,p) > 0) .* (vertexwise(:,v)>0));
        pctoverlap(p) = nnz(parcelwise(parcelinds,p)==vertexwise(parcelinds,v)) / nnz(parcelinds);
        for parcelnum = 1:length(parcelIDs);
            parcelinds = logical(parcels==parcelIDs(parcelnum));
            if (mean(parcelwise(parcelinds,p)) > 0) && (any(vertexwise(parcelinds,v)>0))
                parcelmatchpct = nnz(parcelwise(parcelinds,p)==vertexwise(parcelinds,v)) / nnz(parcelinds .* (vertexwise(:,v)>0));
            end
            overlap(parcelinds,p) = parcelmatchpct;
        end
        overlap(logical(parcelwise(:,p) < 1),p) = -1;
    end
    
%     for p = 1:size(parcelwise,2)
%         parcelmatchpcts = zeros(length(parcelIDs),1);
%         for parcelnum = 1:length(parcelIDs);
%             parcelinds = logical(parcels==parcelIDs(parcelnum));
%             if mean(parcelwise(parcelinds,p)) > 0;
%                 parcelmatchpcts(parcelnum) = nnz(parcelwise(parcelinds,p)==vertexwise(parcelinds,v)) / nnz(parcelinds);
%             end
%         end
%         %parcelmatchpcts(parcelmatchpcts==0) = [];
%             
%         pctoverlap(p) = mean(parcelmatchpcts); 
%     end
    
    
%     [bestmaxpct(v) maxi] = max(pctoverlap);
%     output(:,v) = parcelwise(:,maxi);
%     disp(['Vertexwise thresh ' num2str(vertthreshs(v)) ' best overlap: parcel thresh ' num2str(realparcelthreshs(maxi)) ', overlap = ' num2str(bestmaxpct(v))])
%     
%     threshsoutsidedist = logical(abs(realparcelthreshs-vertthreshs(v)) > maxthreshdist);
%     pctoverlap(threshsoutsidedist) = 0;
%     [bestnearbymaxpct(v) maxi] = max(pctoverlap);
%     output2(:,v) = parcelwise(:,maxi);
%     disp(['Vertexwise thresh ' num2str(vertthreshs(v)) ' best overlap within ' num2str(maxthreshdist) ': parcel thresh ' num2str(realparcelthreshs(maxi)) ', overlap = ' num2str(bestnearbymaxpct(v))])
%     
    %parcelwise(:,maxi) = 1000;
    matchingthreshnum = find(single(realparcelthreshs)==single(vertthreshs(v)));
    
    if isempty(matchingthreshnum)
        disp(['no matching thresh; using smallest'])
        matchingthreshnum = 1;
    end
    vert_pctoverlap(v) = pctoverlap(matchingthreshnum);
    disp(['Vertexwise thresh ' num2str(vertthreshs(v)) ' same thresh: overlap = ' num2str(vert_pctoverlap(v))])
    if ~isempty(matchingthreshnum)
        output3(:,v) = parcelwise(:,matchingthreshnum);
        overlapoutput(:,v) = overlap(:,matchingthreshnum);
    end



end
 templatefile = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';
% cifti_write_wHDR(output,templatefile,'120_Consensus_LR_allthresh_matchingparcelthreshs')
% cifti_write_wHDR(output2,templatefile,'120_Consensus_LR_allthresh_nearbyparcelthreshs')
%cifti_write_wHDR(output3,templatefile,'120_Consensus_LR_allthresh_sameparcelthreshs')
cifti_write_wHDR(overlapoutput,templatefile,'120_Consensus_LR_allthresh_overlapwithsameparcelthreshs')

% save(gifti(single(output)),'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_vertexwise_infomap/geodist/120_Consensus_L_allthresh_matchingparcelthreshs.func.gii')
% save(gifti(single(output2)),'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_vertexwise_infomap/geodist/120_Consensus_L_allthresh_nearbyparcelthreshs.func.gii')

%save(gifti(single(output)),'/data/cn4/evan/RestingState/FC_Mapping_120/120_R_vertexwise_infomap/eucdist/120_Consensus_R_allthresh_matchingparcelthreshs.func.gii')
%save(gifti(single(output2)),'/data/cn4/evan/RestingState/FC_Mapping_120/120_R_vertexwise_infomap/eucdist/120_Consensus_R_allthresh_nearbyparcelthreshs.func.gii')

    

