subname = 'Poldrome';%'MSC02';
cifti_file = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_LR_timeseries.dtseries.nii';%'/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_LR_timeseries.dtseries.nii';
tmaskfile = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_total_tmask.txt';%'/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_total_tmask.txt';
tmask = load(tmaskfile);
surfaceareafiles = {'poldrack.L.midthickness.32k_fs_LR_surfaceareas.func.gii','poldrack.L.midthickness.32k_fs_LR_surfaceareas.func.gii'};

xdistance = 30;

minclustersizemm = 20;

%kden_thresh = .1;
kden_threshs = [.02 : .01 : .1];
r_thresh = .2;
r_thresh_templates = 0;


%-----------------------------

load('/data/cn4/evan/RestingState/Ind_variability/Templates_consensus.mat')
templates{1}(:,11:12) = [];
IDs{1}(11:12) = [];
templates{1}(:,end-1:end) = [];
IDs{1}(end-1:end) = [];

ncortverts = 29696 + 29716;

load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
distances = distances(1:ncortverts,1:ncortverts);

threshdistance = distances > xdistance;

clear distances

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


for kden_thresh = kden_threshs;
disp(kden_thresh)

values_sorted = sort(templates{1}(:),'descend');
threshval = values_sorted(round(numel(templates{1}) .* kden_thresh));
ThreshTemplates = templates{thresh} >= threshval;


subdata = cifti_read(cifti_file);
subdata = subdata(:,logical(tmask));

correlmaps = paircorr_mod(subdata(1:ncortverts,:)');
correlmaps(isnan(correlmaps)) = 0;
correlmaps = FisherTransform(correlmaps);
[i r_sub_thresh kden] = matrix_thresholder_faster2(correlmaps,kden_thresh,'kden');

thissub_networkconnections = zeros(ncortverts,1);

for i = 1:ncortverts
    inds = threshdistance(:,i);
    
    Threshvertmap = correlmaps(i,:) >= r_sub_thresh;

        for templatenum = 1:size(ThreshTemplates,2);
            dice_coeffs(templatenum) = nnz(ThreshTemplates(inds,templatenum) .* Threshvertmap(inds)') ./ nnz(ThreshTemplates(inds,templatenum) + Threshvertmap(inds)');
        end
    
    [maxdice maxi] = max(dice_coeffs);
    
    thissub_networkconnections(i) = IDs{thresh}(maxi);
    
end
clear dice_coeffs



temp = zeros(ncortverts,1);
for ID = IDs{thresh}(:)'
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
cifti_write_wHDR(out,[],['Templatematch_dice_' subname '_kden_' num2str(kden_thresh)]);
end













