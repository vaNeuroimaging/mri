
MSCnames = {'MSC01','MSC02','MSC03','MSC04'};


load('/data/cn4/evan/RestingState/Ind_variability/Templates_consensus.mat')
templates{1}(:,11:12) = [];
IDs{1}(11:12) = [];
%templates{1}(:,11) = [];
%IDs{1}(11) = [];
templates{1}(:,end-1:end) = [];
IDs{1}(end-1:end) = [];

kden_thresh = .1;

values_sorted = sort(templates{1}(:),'descend');
threshval = values_sorted(round(numel(templates{1}) .* kden_thresh));

ncortverts = 29696 + 29716;

xdistance = 30;

minclustersizemm = 50;


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

templatethresh = 1;

%%

load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
distances = distances(1:ncortverts,1:ncortverts);
threshdistance = distances > xdistance;

clear distances
%ThreshTemplates = zeros(size(templates{templatethresh}));
networkconnections = zeros(ncortverts,length(MSCnames));


% for templatenum = 1:size(templates{templatethresh},2);
%     sortvals = sort(templates{templatethresh}(:,templatenum),'descend');
%     threshval = sortvals(round(numel(sortvals)*kden_thresh));
%     ThreshTemplates(:,templatenum) = templates{templatethresh}(:,templatenum) >= threshval;
% end

ThreshTemplates = templates{templatethresh} >= threshval;

for s = 3%1:length(MSCnames)
    disp(['Subject ' num2str(s)])
    MSCname = MSCnames{s};
    tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
    [subjects tmasks] = textread(tmaskfile,'%s %s');
    for i = 1:length(subjects)
        if i == 1
            subdata = cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
            tmask = load(tmasks{i});
        else
            subdata = [subdata cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'])];
            tmask = [tmask; load(tmasks{i})];
        end
    end
    
    subdata = subdata(:,logical(tmask));
    
    
    hems = {'L','R'};
    for hemnum = 1:length(hems)
        surfaceareafiles{hemnum} = ['/data/nil-bluearc/GMT/Laumann/MSC/freesurfer/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/' MSCname '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
        if ~exist(surfaceareafiles{hemnum})
            surffile = [surfaceareafiles{hemnum}(1:end-22) '.surf.gii'];
            %system(['wb_command -surface-vertex-areas ' surffile ' ' surfaceareafiles{hemnum}])
            slashlocs = strfind(surfaceareafiles{hemnum},'/');
            system(['wb_command -surface-vertex-areas ' surffile ' ' surfaceareafiles{hemnum}(slashlocs(end)+1 : end)])
            input_surfaceareafiles{hemnum} = surfaceareafiles{hemnum}(slashlocs(end)+1 : end);
        end
    end
    
    
    
    correlmaps = paircorr_mod(subdata(1:ncortverts,:)');
    correlmaps(isnan(correlmaps)) = 0;
    correlmaps = FisherTransform(correlmaps);
    
    thissub_networkconnections = zeros(ncortverts,1);
    
    for i = 1:ncortverts

        inds = threshdistance(:,i);
        
        for templatenum = 1:size(templates{templatethresh},2);
            
            meanvals(templatenum) = mean(correlmaps(logical(inds .* ThreshTemplates(1:ncortverts,templatenum)),i));
            
        end
        
        [maxdice maxi] = max(meanvals);
        
        thissub_networkconnections(i) = IDs{templatethresh}(maxi);
        
    end
    clear dice_coeffs
    
    temp = zeros(ncortverts,1);
    for ID = IDs{templatethresh}(:)'
        %outputcifti = metric_cluster_cifti_surfacearea(thissub_networkconnections,ID-.5,ID+.5,minclustersizemm,surfaceareafiles);
        outputcifti = metric_cluster_cifti_surfacearea(thissub_networkconnections,ID-.5,ID+.5,minclustersizemm,input_surfaceareafiles);
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
    
    
    clear correlmaps subdata outputcifti thissub_networkconnections
end
out = zeros(66697,size(networkconnections,2));
out(1:size(networkconnections,1),:) = networkconnections;
cifti_write_wHDR(out,[],['Templatematch_meanvals_bysubject_' num2str(kden_thresh) 'kden_sub3']);


%%

% subsizes_totest = [20 30 40 50 60 70];
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





