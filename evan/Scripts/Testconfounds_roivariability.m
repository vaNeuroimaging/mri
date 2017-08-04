clustersizefileL = '/data/cn4/evan/RestingState/Ind_variability/120/Variability_L_consensus_dice_clustersize.func.gii';
clustersizefileR = '/data/cn4/evan/RestingState/Ind_variability/120/Variability_R_consensus_dice_clustersize.func.gii';

mask{1} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); mask{1} = ~mask{1}.cdata;
mask{2} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); mask{2} = ~mask{2}.cdata;

% clustersizeL = gifti(clustersizefileL);
% clustersizeL_thresh = clustersizeL.cdata(logical(maskL),2) > clustersizethresh;
% clustersizeR = gifti(clustersizefileR);
% clustersizeR_thresh = clustersizeR.cdata(logical(maskR),2) > clustersizethresh;
%
% variable_regions = zeros(66697,1);
% variable_regions(1:(nnz(maskL)+nnz(maskR))) = [clustersizeL_thresh;clustersizeR_thresh];
% cifti_write_wHDR(variable_regions,[],['Variable_regions_' num2str(clustersizethresh)])
%

%gifti_to_cifti('Variabile_regions_L_consensus_dice_distance10.func.gii','Variabile_regions_R_consensus_dice_distance10.func.gii','Variabile_regions_LR_consensus_dice_distance10')

all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance10.dtseries.nii');

%networkIDtotest = 8;
%all_variable_regions(all_variable_regions~=networkIDtotest) = 0;

variable_regions = metric_cluster_cifti(all_variable_regions,.5,100,0); 



tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

all_diffcurv = zeros(size(all_variable_regions,1),length(subjects));
all_diffsulc = zeros(size(all_variable_regions,1),length(subjects));

hems = {'L','R'};

for s = 1:length(tmasks)
tmask = load(tmasks{s});
ntimepoints(s) = nnz(tmask);

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    subcurv = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' hem '.curvature.32k_fs_LR.shape.gii']);
    subsulc = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' hem '.sulc.32k_fs_LR.shape.gii']);
    groupcurv = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.curvature.32k_fs_LR.shape.gii']);
    groupsulc = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sulc.32k_fs_LR.shape.gii']);
    diffcurv = subcurv.cdata - groupcurv.cdata;
    diffsulc = subsulc.cdata - groupsulc.cdata;
    
    inds = [1:nnz(mask{hemnum})] + (nnz(mask{1}) .* (hemnum-1));

    all_diffcurv(inds,s) = diffcurv(logical(mask{hemnum}));
    all_diffsulc(inds,s) = diffsulc(logical(mask{hemnum}));
end

end

%variable_regions(:,1) = [];

%%
% load('/data/cn4/evan/RestingState/Ind_variability/120/Templates_consensus.mat')
% templates{1}(:,11:12) = [];
% IDs{1}(11:12) = [];
% templates{1}(:,end-1:end) = [];
% IDs{1}(end-1:end) = [];
% 
% templatefile = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';
% 
% kden_thresh = .1;
% 
% ncortverts = nnz(maskL)+nnz(maskR);
% 
% load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
% distances = distances(1:ncortverts,1:ncortverts);
% 
% tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_TMASKLIST.txt';
% [subjects tmasks] = textread(tmaskfile,'%s %s');
% 
% xdistance = 20;
% 
% minvertclustersize = 10;
% 
% threshdistance = distances > xdistance;
% 
% clear distances
% 
% xdist_inds = zeros(size(variable_regions,2),size(threshdistance,1));
% for i = 1:size(variable_regions,2)
%     xdist_inds(i,:) = logical(all(threshdistance(logical(variable_regions(:,i)),:),1));
% end
% 
% xdist_inds = logical(xdist_inds);
% 
% clear threshdistance
% 
% thresh = 1;
% ThreshTemplates = zeros(size(templates{thresh}));
% 
% 
% 
% for templatenum = 1:size(templates{thresh},2);
%     sortvals = sort(templates{thresh}(:,templatenum),'descend');
%     threshval = sortvals(round(numel(sortvals)*kden_thresh));
%     ThreshTemplates(:,templatenum) = templates{thresh}(:,templatenum) >= threshval;
% end
% 
% networkconnections = zeros(size(variable_regions,2),length(subjects));
% 
% for s = 1:length(subjects)
%     disp(['Subject ' num2str(s)])
%     subject = subjects{s};
%     cifti_file = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subject '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
%     subdata = cifti_read(cifti_file);
%     tmask = load(tmasks{s});
%     subdata = subdata(:,logical(tmask));
%     
%     
%     for i = 1:size(variable_regions,2)
%         timecourse = mean(subdata(logical(variable_regions(:,i)),:),1);
%         correlmap = paircorr_mod(timecourse',subdata(1:ncortverts,:)');
%         correlmap(isnan(correlmap)) = 0;
%         correlmap = FisherTransform(correlmap);
%         
%         sortvals = sort(correlmap,'descend');
%         threshval = sortvals(round(numel(correlmap)*kden_thresh));
%         Threshvertmap = correlmap >= threshval;
%         
%         %inds = logical(all(threshdistance(logical(variable_regions(:,i)),:),1));
%         
%         for templatenum = 1:size(templates{thresh},2);
%             dice_coeffs(templatenum) = nnz(ThreshTemplates(xdist_inds(templatenum,:),templatenum) == Threshvertmap(xdist_inds(templatenum,:))') ./ nnz(xdist_inds(templatenum,:));
%         end
%         
%         [maxdice maxi] = max(dice_coeffs);
%         
%         networkconnections(i,s) = IDs{thresh}(maxi);
%         
%     end
%     clear dice_coeffs
%     
%     clear correlmaps subdata outputcifti thissub_networkconnections
% end

%%

surfaceareas = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject.dtseries.nii');

networkconnections = zeros(size(variable_regions,2),size(networkconnection_bysub,2));
networkconnections_SAincluster = zeros(size(variable_regions,2),size(networkconnection_bysub,2));
curvincluster = zeros(size(variable_regions,2),size(networkconnection_bysub,2));
sulcincluster = zeros(size(variable_regions,2),size(networkconnection_bysub,2));

for i = 1:size(variable_regions,2)
    
    targetID = mode(variable_regions(logical(variable_regions(:,i)),i));
    
    for s = 1:size(networkconnection_bysub,2)
        
        %networkconnections(i,s) = mode(networkconnection_bysub(logical(variable_regions(:,i)),s));
        
        if sum(surfaceareas(logical((networkconnection_bysub(:,s)==targetID) .* (variable_regions(:,i))))) >= (.25 * sum(surfaceareas(logical(variable_regions(:,i)))))
            %nnz(networkconnection_bysub(logical(variable_regions(:,i)),s)==targetID) >= (.05 * nnz(variable_regions(:,i)))
            %networkconnections(i,s) = targetID;
            networkconnections(i,s) = 1;
        else
            %networkconnections(i,s) = mode(networkconnection_bysub(logical(variable_regions(:,i)),s));
            networkconnections(i,s) = 0;
        end
        
        networkconnections_SAincluster(i,s) = sum(surfaceareas(logical((networkconnection_bysub(:,s)==targetID) .* (variable_regions(:,i)))));
        curvincluster(i,s) = sum(abs(all_diffcurv(logical(variable_regions(:,i)),s)));
        sulcincluster(i,s) = sum(abs(all_diffsulc(logical(variable_regions(:,i)),s)));
        
    end
    
end
%%
disp('Surface Area with alternate connection corsscorrelations, controlling for # timepoints')
[R,P]=partialcorr(networkconnections_SAincluster',ntimepoints');
disp(R)
disp(P .* nnz(triu(P,1)))

curvR = zeros(1,size(variable_regions,2));
curvP = zeros(1,size(variable_regions,2));
sulcR = zeros(1,size(variable_regions,2));
sulcP = zeros(1,size(variable_regions,2));
for i = 1:size(variable_regions,2)
    [curvR(i) curvP(i)] = corr(networkconnections_SAincluster(i,:)',curvincluster(i,:)');
    [sulcR(i) sulcP(i)] = corr(networkconnections_SAincluster(i,:)',sulcincluster(i,:)');
end

disp('Surface Area with alternate connection vs curvature')
disp(curvR)
disp(curvP .* numel(curvP))

disp('Surface Area with alternate connection vs sulcal depth')
disp(sulcR)
disp(sulcP .* numel(sulcP))