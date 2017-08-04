sizethresh = 10;
clustersizefileL = '/data/cn4/evan/RestingState/Ind_variability/120/Variability_L_consensus_dice_clustersize.func.gii';
clustersizefileR = '/data/cn4/evan/RestingState/Ind_variability/120/Variability_R_consensus_dice_clustersize.func.gii';

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;

% clustersizeL = gifti(clustersizefileL);
% clustersizeL_thresh = clustersizeL.cdata(logical(maskL),2) > clustersizethresh;
% clustersizeR = gifti(clustersizefileR);
% clustersizeR_thresh = clustersizeR.cdata(logical(maskR),2) > clustersizethresh;
% 
% variable_regions = zeros(66697,1);
% variable_regions(1:(nnz(maskL)+nnz(maskR))) = [clustersizeL_thresh;clustersizeR_thresh];
% cifti_write_wHDR(variable_regions,[],['Variable_regions_' num2str(clustersizethresh)])
% 


variable_regions = cifti_read('Variable_regions_dist_LR_VA.dtseries.nii');
variable_regions = sum(variable_regions,2);
variable_inds = find(variable_regions);

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
% xdistance = 30;
% 
% minvertclustersize = 10;
% 
% threshdistance = distances > xdistance;
% 
% clear distances
% 
% thresh = 1;
%     ThreshTemplates = zeros(size(templates{thresh}));
%     
%     
%     
%     for templatenum = 1:size(templates{thresh},2);
%         sortvals = sort(templates{thresh}(:,templatenum),'descend');
%         threshval = sortvals(round(numel(sortvals)*kden_thresh));
%         ThreshTemplates(:,templatenum) = templates{thresh}(:,templatenum) >= threshval;
%     end
%     
%     networkconnections = zeros(length(variable_inds),length(subjects));
%     
%     for s = 1:length(subjects)
%         disp(['Subject ' num2str(s)])
%         subject = subjects{s};
%         cifti_file = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subject '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
%         subdata = cifti_read(cifti_file);
%         tmask = load(tmasks{s});
%         subdata = subdata(:,logical(tmask));
%         
%         correlmaps = paircorr_mod(subdata(variable_inds,:)',subdata(1:ncortverts,:)');
%         correlmaps(isnan(correlmaps)) = 0;
%         correlmaps = FisherTransform(correlmaps);
%         
%         
%         
%         for inum = 1:length(variable_inds)
%             i = variable_inds(inum);
%             inds = threshdistance(:,i);
%             
%             sortvals = sort(correlmaps(inum,:),'descend');
%             threshval = sortvals(round(numel(correlmaps(inum,:))*kden_thresh));
%             Threshvertmap = correlmaps(inum,:) >= threshval;
%                 
%             for templatenum = 1:size(templates{thresh},2);
%                 dice_coeffs(templatenum) = nnz(ThreshTemplates(inds,templatenum) == Threshvertmap(inds)') ./ nnz(inds);
%             end
%             
%             [maxdice maxi] = max(dice_coeffs);
%             
%             networkconnections(inum,s) = IDs{thresh}(maxi);
% 
%             
%         end
%         clear dice_coeffs
%         
%         clear correlmaps subdata outputcifti thissub_networkconnections
%     end
    
    %%
    
    
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject.dtseries.nii');

networkconnections = zeros(size(variable_regions,2),size(networkconnection_bysub,2));

for inum = 1:length(variable_inds)
    for s = 1:size(networkconnection_bysub,2)
        
        networkconnections(inum,s) = networkconnection_bysub(variable_inds(inum),s);
        
    end
    
end

    
Y = pdist(networkconnections','hamming');           % 'pdist' converts the square adjacency matrix to a
%  1 x n matrix so that the function linkage can construct the tree

clustering = linkage(squareform(Y), 'average','hamming');      % 'linkage' computes the data to construct the tree
% 'average' refers to the UPGMA algorithm

clusters = cluster(clustering, 'MaxClust', [1:80]);

for s = 1:size(networkconnections,2)
    for s2 = 1:size(networkconnections,2)
        similaritymat(s,s2) = nnz(networkconnections(:,s)==networkconnections(:,s2));
    end
end


for numclust = 1:size(clusters,2)
    Qvals(numclust) = M_calc_modularity(clusters(:,numclust),similaritymat);
end
[maxQval maxQind] = max(Qvals);

cut = clustering(end-maxQind,3);


colormap = [1 0 0;0 0 1;1 1 0; 0 .8 0; 1 .6 1;1 .5 0;1 .7 .4;0 .6 .6;.6 .2 1];
    %0 0 0;.2 1 1;;0 .2 .4;.2 1 .2;0 0 1;1 1 .8;0 .4 0;.25 .25 .25];

[H,T,perm] = dendrogram_evan(clustering, colormap,sizethresh,0, 'orientation','left', 'colorthreshold', cut+.00001 ,'ORIENTATION','left');
        
%[H,T,perm] = dendrogram(clustering, 0, 'orientation','left', 'colorthreshold', cut,'ORIENTATION','left');

cophenetic_r = cophenet(clustering, Y);
disp(cophenetic_r)

%%
clear clustersizes IDs

clusters = cluster(clustering,'Cutoff',cut+.00001,'Criterion','distance');
clusternums = unique(clusters);
IDs = zeros(size(networkconnections,1),0);
clustersizes = [];
for i = 1:length(clusternums)
    if nnz(clusters==clusternums(i)) >= sizethresh
        clustersizes(end+1) = nnz(clusters==clusternums(i));
        IDs(:,end+1) = mode(networkconnections(:,clusters==clusternums(i)),2);
    end
end

out = zeros(66697,size(IDs,2));
for i = 1:size(IDs,1)
    for j = 1:size(IDs,2)
        out(variable_inds(i),j) = IDs(i,j);
    end
end
clustersizes
pause(.5)
cifti_write_wHDR(out,[],'Variable_regions_dist_LR_clustered_VA')
    
    

