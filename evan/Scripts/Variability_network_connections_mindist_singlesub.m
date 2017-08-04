subname = 'MSC02';
cifti_file = '/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_LR_timeseries.dtseries.nii';
tmaskfile = '/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_total_tmask.txt';
tmask = load(tmaskfile);
surfaceareafiles = {'MSC02.L.midthickness.32k_fs_LR_surfaceareas_func.gii','MSC02.R.midthickness.32k_fs_LR_surfaceareas_func.gii'};

xdistance = 30;

minclustersizemm = 50;

%kden_thresh = .1;
%kden_threshs = [.005 : .005 : .045];
kden_threshs = .1;
r_thresh = .2;
r_thresh_templates = 0;


%-----------------------------

ncortverts = 29696 + 29716;

groupnetworksfile = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR.dtseries.nii';
groupnetworks = cifti_read(groupnetworksfile); groupnetworks = groupnetworks(1:ncortverts);
groupnetworkIDs = unique(groupnetworks); groupnetworkIDs(groupnetworkIDs==0) = [];

% load('/data/cn4/evan/RestingState/Ind_variability/120/Templates_consensus.mat')
% templates{1}(:,11:12) = [];
% IDs{1}(11:12) = [];
% templates{1}(:,end-1:end) = [];
% IDs{1}(end-1:end) = [];



%  load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
%  distances = distances(1:ncortverts,1:ncortverts);
 
 %threshdistance = distances > xdistance;


bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

medial_wall{1} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);
medial_wall{1} = medial_wall{1}.cdata;
medial_wall{2} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);
medial_wall{2} = medial_wall{2}.cdata;

thresh = 1;


% subdata = cifti_read(cifti_file);
% subdata = subdata(:,logical(tmask));
% 
% correlmaps = paircorr_mod(subdata(1:ncortverts,:)');
% correlmaps(isnan(correlmaps)) = 0;
% correlmaps = FisherTransform(correlmaps);



for groupnetworknum = 1:length(groupnetworkIDs)
    groupnetworkinds{groupnetworknum} = find(groupnetworks==groupnetworkIDs(groupnetworknum));
    groupnetworkdists{groupnetworknum} = distances(:,groupnetworkinds{groupnetworknum});
end

for kden_thresh = kden_threshs;


% ThreshTemplates = zeros(size(templates{thresh}));
% 
% 
% for templatenum = 1:size(templates{thresh},2);
%     if r_thresh_templates
%         ThreshTemplates(:,templatenum) = templates{thresh}(:,templatenum) >= r_thresh;
%     else
%         sortvals = sort(templates{thresh}(:,templatenum),'descend');
%         threshval = sortvals(round(numel(sortvals)*kden_thresh));
%         ThreshTemplates(:,templatenum) = templates{thresh}(:,templatenum) >= threshval;
%     end
% end




thissub_networkconnections = zeros(ncortverts,1);
%matlabpool open 4
for i = 1:ncortverts
    disp(num2str(i))
    distinds = find(distances(:,i) > xdistance);
    
    sortvals = sort(correlmaps(i,:),'descend');
    threshval = sortvals(round(numel(correlmaps(i,:))*kden_thresh));
    Threshvertinds = intersect(find(correlmaps(i,:)' >= threshval),distinds);
    
%     node_groupnetworks = groupnetworks_separated .* repmat(distinds,1,size(groupnetworks_separated,2));
%     
%     1;
    mean_mindists = zeros(length(groupnetworkIDs),1);
    for groupnetworknum = 1:length(groupnetworkIDs)
        
        thisnetworkpattern = groupnetworkdists{groupnetworknum};
        
        [ign groupinds_touse ign2] = intersect(groupnetworkinds{groupnetworknum},distinds);
        thesedistances = thisnetworkpattern(Threshvertinds,groupinds_touse);
        
        groupmindists = min(thesedistances,[],2);
        groupmindists(groupmindists>50) = [];
        
        submindists = min(thesedistances,[],1);
        submindists(submindists > 50) = [];
        
        mean_mindists(groupnetworknum) = mean([submindists(:); groupmindists(:)]);
        
        
    end
    
    [mindist mini] = min(mean_mindists);
    
    thissub_networkconnections(i) = groupnetworkIDs(mini);
    
end
%matlabpool close
%profile viewer


temp = zeros(ncortverts,1);
for ID = groupnetworkIDs(:)'
    outputcifti = metric_cluster_cifti_surfacearea(thissub_networkconnections,ID-.5,ID+.5,minclustersizemm,surfaceareafiles);
    outputcifti(logical(outputcifti)) = ID;
    temp = temp + sum(outputcifti,2);
end

thissub_networkconnections = temp;

for hem = 1:2
    gifti_networkconnections = zeros(32492,1);
    gifti_networkconnections(medial_wall{hem}==0) = thissub_networkconnections((1:nnz(medial_wall{hem}==0)) + (nnz(medial_wall{1}==0) * (hem-1)));
    
    blankinds = find((gifti_networkconnections==0) .* (medial_wall{hem}==0));
    while ~isempty(blankinds)
        temp = gifti_networkconnections;
        for ind = blankinds(:)'
            indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
            neighvals = gifti_networkconnections(indneighs); neighvals(neighvals==0) = [];
            if ~isempty(neighvals)
                temp(ind) = mode(neighvals);
            end
        end
        gifti_networkconnections = temp;
        blankinds = find((gifti_networkconnections==0) .* (medial_wall{hem}==0));
    end
    
    thissub_networkconnections((1:nnz(medial_wall{hem}==0)) + (nnz(medial_wall{1}==0) * (hem-1))) = gifti_networkconnections(medial_wall{hem}==0);
    
end

out = zeros(66697,1);
out(1:size(thissub_networkconnections,1)) = thissub_networkconnections;
cifti_write_wHDR(out,[],['Templatematch_mindist_' subname '_kden_' num2str(kden_thresh)]);
end













