% subname = 'Poldrome';%'MSC02';
% cifti_file = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_LR_timeseries.dtseries.nii';%'/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_LR_timeseries.dtseries.nii';
% tmaskfile = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_total_tmask.txt';%'/net/nil-bluearc/GMT/Laumann/MSC/MSC02/Functionals/FCPROCESS_SCRUBBED/cifti_correlation_normalwall/10ses_concat/allsubs_total_tmask.txt';
% tmask = load(tmaskfile);
% surfaceareafiles = {'poldrack.L.midthickness.32k_fs_LR_surfaceareas.func.gii','poldrack.L.midthickness.32k_fs_LR_surfaceareas.func.gii'};


tmaskfile = '/data/cn4/evan/Nelson/VISN17_tmasklist.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');


subnums = [1:length(subjects)];%[121:228];

subjects = subjects(subnums);
tmasks = tmasks(subnums);
%ciftifiles = ciftifiles(subnums);

xdistance = 30;

minclustersizemm = 30;

%kden_thresh = .1;
kden_threshs = .05;%[.02 : .01 : .1];
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
    %ThreshTemplates = zeros(size(templates{thresh}));
    
    
    % for templatenum = 1:size(templates{thresh},2);
    %     if r_thresh_templates
    %         ThreshTemplates(:,templatenum) = templates{thresh}(:,templatenum) >= r_thresh;
    %     else
    %         sortvals = sort(templates{thresh}(:,templatenum),'descend');
    %         threshval = sortvals(round(numel(sortvals)*kden_thresh));
    %         ThreshTemplates(:,templatenum) = templates{thresh}(:,templatenum) >= threshval;
    %     end
    % end
    
    
    for s = 1:length(subjects)
        disp(['Subject ' num2str(s)])
        subject = subjects{s};
        %cifti_file = ciftifiles{s};
        cifti_file = ['/data/cn4/evan/Nelson/cifti_normalwall_timeseries/' subject '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
        
        hems = {'L','R'};
        for hemnum = 1:length(hems)
            surfaceareafiles{hemnum} = ['/data/cn4/evan/Nelson/FREESURFER/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
            if ~exist(surfaceareafiles{hemnum})
                surffile = [surfaceareafiles{hemnum}(1:end-22) '.surf.gii'];
                system(['wb_command -surface-vertex-areas ' surffile ' ' surfaceareafiles{hemnum}])
            end
        end
        
        subdata = cifti_read(cifti_file);
        tmask = load(tmasks{s});
        subdata = subdata(:,logical(tmask));
        
        correlmaps = paircorr_mod(subdata(1:ncortverts,:)');
        correlmaps(isnan(correlmaps)) = 0;
        correlmaps = FisherTransform(correlmaps);
        
        [i r_sub_thresh kden] = matrix_thresholder_faster2(correlmaps,kden_thresh,'kden');
        
        thissub_networkconnections = zeros(ncortverts,1);
        
        for i = 1:ncortverts
            inds = threshdistance(:,i);
            
            %     sortvals = sort(correlmaps(i,:),'descend');
            %     threshval = sortvals(round(numel(correlmaps(i,:))*kden_thresh));
            %     Threshvertmap = correlmaps(i,:) >= threshval;
            %
            %     for templatenum = 1:size(templates{thresh},2);
            %         dice_coeffs(templatenum) = nnz(ThreshTemplates(inds,templatenum) == Threshvertmap(inds)') ./ nnz(inds);
            %     end
            
            Threshvertmap = correlmaps(i,:) >= r_sub_thresh;
            if nnz(Threshvertmap)<10
                sortedvals = sort(correlmaps(i,:),'descend');
                r_temp_thresh = sortedvals(round(kden_thresh*size(correlmaps,1)));
                Threshvertmap = correlmaps(i,:) >= r_temp_thresh;
            end
            
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
        networkconnections(:,s) = thissub_networkconnections;
    end
    
    out = zeros(66697,size(networkconnections,2));
    out(1:size(networkconnections,1),:) = networkconnections;
    cifti_write_wHDR(out,[],['Templatematch_dice_bysubject_kden' num2str(kden_thresh)]);
    
end

%%
% 
% templatethresh = 1;
% 
% subsizes_totest = round([.1*length(subjects) .15*length(subjects) .2*length(subjects) .25*length(subjects) .3*length(subjects)  .35*length(subjects) .4*length(subjects)]);
% 
% maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = maskL.cdata;
% maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = maskR.cdata;
% 
% IDorder = zeros(size(networkconnections,1),length(IDs{templatethresh}));
% IDsize = zeros(size(networkconnections,1),length(IDs{templatethresh}));
% 
% for i = 1:ncortverts
%     
%     thisvert_IDs = networkconnections(i,:);
%     
%     for IDnum = 1:length(IDs{templatethresh})
%         IDcounts(IDnum) = nnz(thisvert_IDs==IDs{templatethresh}(IDnum));
%     end
%     
%     [IDcounts_sorted sorti] = sort(IDcounts,'descend');
%     clear IDcounts
%     
%     IDorder(i,1:length(sorti)) = IDs{templatethresh}(sorti);
%     IDsize(i,1:length(sorti)) = IDcounts_sorted;
%     
% end
% 
% ID_order_L = zeros(size(maskL,1),length(IDs{templatethresh}));
% ID_order_L(maskL==0,:) = IDorder(1:nnz(maskL==0),:);
% 
% ID_size_L = zeros(size(maskL,1),length(IDs{templatethresh}));
% ID_size_L(maskL==0,:) = IDsize(1:nnz(maskL==0),:);
% 
% ID_order_R = zeros(size(maskR,1),length(IDs{templatethresh}));
% ID_order_R(maskR==0,:) = IDorder(nnz(maskL==0)+1 : nnz(maskR==0) + nnz(maskL==0),:);
% 
% ID_size_R = zeros(size(maskR,1),length(IDs{templatethresh}));
% ID_size_R(maskR==0,:) = IDsize(nnz(maskL==0)+1 : nnz(maskR==0) + nnz(maskL==0),:);
% 
% %     save(gifti(single(ID_order_L)),['Variability_L_thresh' num2str(thresh) '.func.gii'])
% %     save(gifti(single(ID_size_L)),['Variability_L_thresh' num2str(thresh) '_clustersize.func.gii'])
% %     Cluster_and_combine(['Variability_L_thresh' num2str(thresh) '.func.gii'],'L')
% %
% %     save(gifti(single(ID_order_R)),['Variability_R_thresh' num2str(thresh) '.func.gii'])
% %     save(gifti(single(ID_size_R)),['Variability_R_thresh' num2str(thresh) '_clustersize.func.gii'])
% %     Cluster_and_combine(['Variability_R_thresh' num2str(thresh) '.func.gii'],'R')
% 
% save(gifti(single(ID_order_L)),['Variability_L_consensus_dice.func.gii'])
% save(gifti(single(ID_size_L)),['Variability_L_consensus_dice_clustersize.func.gii'])
% 
% 
% save(gifti(single(ID_order_R)),['Variability_R_consensus_dice.func.gii'])
% save(gifti(single(ID_size_R)),['Variability_R_consensus_dice_clustersize.func.gii'])
% 
% for sizethresh = subsizes_totest
%     Cluster_and_combine(['Variability_L_consensus_dice.func.gii'],'L',sizethresh)
%     Cluster_and_combine(['Variability_R_consensus_dice.func.gii'],'R',sizethresh)
% end
% 
% 










