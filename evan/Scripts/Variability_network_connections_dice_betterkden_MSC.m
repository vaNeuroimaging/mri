%cifti_file = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_LR_timeseries.dtseries.nii';%'/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_LR_timeseries.dtseries.nii';
%tmaskfile = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_total_tmask.txt';%'/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_total_tmask.txt';
%tmask = load(tmaskfile);
MSCname = 'MSC01';

surfaceareafiles = {[MSCname '.L.midthickness.32k_fs_LR_surfaceareas.func.gii'],[MSCname '.R.midthickness.32k_fs_LR_surfaceareas.func.gii']};

system(['wb_command -surface-vertex-areas /data/nil-bluearc/GMT/Laumann/MSC/fs5.3/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/' MSCname '.L.midthickness.32k_fs_LR.surf.gii ' surfaceareafiles{1}])
system(['wb_command -surface-vertex-areas /data/nil-bluearc/GMT/Laumann/MSC/fs5.3/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/' MSCname '.R.midthickness.32k_fs_LR.surf.gii ' surfaceareafiles{2}])


medial_wall{1} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);
medial_wall{1} = medial_wall{1}.cdata;
medial_wall{2} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);
medial_wall{2} = medial_wall{2}.cdata;

ncortverts = nnz(~medial_wall{1}) + nnz(~medial_wall{2});


tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
[subjects tmasks] = textread(tmaskfile,'%s %s');
for i = 1:length(subjects)
    if i == 1
        subdata = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
        subdata = subdata.data;
        tmask = load(tmasks{i});
    else
        temp = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
        subdata = [subdata temp.data];
        clear temp
        tmask = [tmask; load(tmasks{i})];
    end
end
subdata = subdata(1:ncortverts,logical(tmask));






xdistance = 30;

minclustersizemm = 50;

kden_threshs = .05;
%kden_threshs = [.02 : .01 : .1];
r_thresh = .2;
r_thresh_templates = 0;


%-----------------------------

load('/data/cn4/evan/RestingState/Ind_variability/Templates_consensus.mat')
templates{1}(:,11:12) = [];
IDs{1}(11:12) = [];
templates{1}(:,end-1:end) = [];
IDs{1}(end-1:end) = [];



bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


threshdistance = zeros(ncortverts);

distances = smartload(['/data/nil-bluearc/GMT/Laumann/MSC/fs5.3/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/geodesic_distance.L.32k_fs_LR.mat']);
distances = distances(~medial_wall{1},~medial_wall{1});
inds = 1:nnz(~medial_wall{1});
threshdistance(inds,inds) = distances;
clear distances

distances = smartload(['/data/nil-bluearc/GMT/Laumann/MSC/fs5.3/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/geodesic_distance.R.32k_fs_LR.mat']);
distances = distances(~medial_wall{2},~medial_wall{2});
inds = (nnz(~medial_wall{1})+1) : (nnz(~medial_wall{1})+nnz(~medial_wall{2}));
threshdistance(inds,inds) = distances;

%load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
%distances = distances(1:ncortverts,1:ncortverts);

%threshdistance = distances > xdistance;

clear distances

threshdistance = threshdistance > xdistance;

thresh = 1;


for kden_thresh = kden_threshs;
disp(kden_thresh)

values_sorted = sort(templates{1}(:),'descend');
threshval = values_sorted(round(numel(templates{1}) .* kden_thresh));
ThreshTemplates = templates{thresh} >= threshval;


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
cifti_write_wHDR(out,[],[MSCname '_Templatematch_dice_kden_' num2str(kden_thresh)]);
end













