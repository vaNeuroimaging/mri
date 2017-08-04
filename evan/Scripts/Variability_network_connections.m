tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_TMASKLIST.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

load('/data/cn4/evan/RestingState/Ind_variability/120/Templates_consensus.mat')
templates{1}(:,11) = [];
IDs{1}(11) = [];

templatefile = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';

minsimilarity = .1;

ncortverts = 29696 + 29716;

load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
distances = distances(1:ncortverts,1:ncortverts);

xdistance = 30;

minvertclustersize = 10;

minsubcluster = 20;

threshdistance = distances > xdistance;

clear distances

for thresh = 1:length(templates)
    FisherTemplates = FisherTransform(templates{thresh});
    networkconnections = zeros(ncortverts,length(subjects));
    
    
    for s = 1:length(subjects)
        disp(['Subject ' num2str(s)])
        subject = subjects{s};
        cifti_file = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subject '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
        subdata = cifti_read(cifti_file);
        tmask = load(tmasks{s});
        subdata = subdata(:,logical(tmask));
        
        correlmaps = paircorr_mod(subdata(1:ncortverts,:)');
        correlmaps(isnan(correlmaps)) = 0;
        correlmaps = FisherTransform(correlmaps)';
        
        thissub_networkconnections = zeros(ncortverts,1);
        
        for i = 1:ncortverts
            inds = threshdistance(:,i);
            
            similarities = paircorr_mod(correlmaps(inds,i),FisherTemplates(inds,:));
            
            [maxcorrel maxi] = max(similarities);
            
            if maxcorrel > minsimilarity
                
                thissub_networkconnections(i) = IDs{thresh}(maxi);
                %networkconnections(i,s) = IDs{thresh}(maxi);
                
            end
        end
        
        networkconnections(:,s) = zeros(ncortverts,1);
        for ID = IDs{thresh}(:)'
            outputcifti = metric_cluster_cifti(thissub_networkconnections,ID-.5,ID+.5,minvertclustersize);
            outputcifti(logical(outputcifti)) = ID;
            networkconnections(:,s) = networkconnections(:,s) + sum(outputcifti,2);
        end
        
        clear correlmaps subdata outputcifti thissub_networkconnections
    end
    
    
    
    
    %%
    
    maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = maskL.cdata;
    maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = maskR.cdata;
    
    IDorder = zeros(size(networkconnections,1),length(IDs{thresh}));
    IDsize = zeros(size(networkconnections,1),length(IDs{thresh}));
    
    for i = 1:ncortverts
        
        thisvert_IDs = networkconnections(i,:);
        
        for IDnum = 1:length(IDs{thresh})
            IDcounts(IDnum) = nnz(thisvert_IDs==IDs{thresh}(IDnum));
        end
        
        [IDcounts_sorted sorti] = sort(IDcounts,'descend');
        clear IDcounts
        
        IDorder(i,1:length(sorti)) = IDs{thresh}(sorti);
        IDsize(i,1:length(sorti)) = IDcounts_sorted;
        
    end
    
    ID_order_L = zeros(size(maskL,1),length(IDs{thresh}));
    ID_order_L(maskL==0,:) = IDorder(1:nnz(maskL==0),:);
    
    ID_size_L = zeros(size(maskL,1),length(IDs{thresh}));
    ID_size_L(maskL==0,:) = IDsize(1:nnz(maskL==0),:);
    
    ID_order_R = zeros(size(maskR,1),length(IDs{thresh}));
    ID_order_R(maskR==0,:) = IDorder(nnz(maskL==0)+1 : nnz(maskR==0) + nnz(maskL==0),:);
    
    ID_size_R = zeros(size(maskR,1),length(IDs{thresh}));
    ID_size_R(maskR==0,:) = IDsize(nnz(maskL==0)+1 : nnz(maskR==0) + nnz(maskL==0),:);
    
%     save(gifti(single(ID_order_L)),['Variability_L_thresh' num2str(thresh) '.func.gii'])
%     save(gifti(single(ID_size_L)),['Variability_L_thresh' num2str(thresh) '_clustersize.func.gii'])
%     Cluster_and_combine(['Variability_L_thresh' num2str(thresh) '.func.gii'],'L')
%     
%     save(gifti(single(ID_order_R)),['Variability_R_thresh' num2str(thresh) '.func.gii'])
%     save(gifti(single(ID_size_R)),['Variability_R_thresh' num2str(thresh) '_clustersize.func.gii'])
%     Cluster_and_combine(['Variability_R_thresh' num2str(thresh) '.func.gii'],'R')

    save(gifti(single(ID_order_L)),['Variability_L_consensus.func.gii'])
    save(gifti(single(ID_size_L)),['Variability_L_consensus_clustersize.func.gii'])
    Cluster_and_combine(['Variability_L_consensus.func.gii'],'L',minsubcluster)
    
    save(gifti(single(ID_order_R)),['Variability_R_consensus.func.gii'])
    save(gifti(single(ID_size_R)),['Variability_R_consensus_clustersize.func.gii'])
    Cluster_and_combine(['Variability_R_consensus.func.gii'],'R',minsubcluster)
    
end









