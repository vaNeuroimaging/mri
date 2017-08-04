subnetworksfile = '';
distthresh = 10;

pctsubthresh = 1/3;

new_cluster = 0;

calc_dices = 0;


subnetworks = ft_read_cifti_mod(subnetworksfile); 
out_template = subnetworks; out_template.data = [];
subnames = subnetworks.mapname;
subnetworks = subnetworks.data;
IDs = unique(subnetworks); IDs(IDs<1) = [];

networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'};

if new_cluster
    
    %surfacearea = ft_read_cifti_mod('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');
    %surfacearea = surfacearea.data;
    
    clusters_subs_IDs_dices_SAs_new = zeros(0,5);
    subnetworks_bycluster_new = zeros(size(subnetworks));
    
    prevstring = [];
    clustercount_new = 0;
    for s = 1:size(subnetworks,2)
        subject = subnames{s};
        string = ['Clustering systems: Subject ' num2str(s) ': ' subject];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        
        SAL = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.midthickness.32k_fs_LR_surfaceareas.func.gii']);
        SAR = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.midthickness.32k_fs_LR_surfaceareas.func.gii']);
        surfacearea = [SAL.cdata; SAR.cdata];
        surfacearea = surfacearea(out_template.brainstructure(1:length(surfacearea))>0);
        
        for ID = IDs(:)'
            outputcifti = metric_cluster_cifti(subnetworks(:,s),ID-.05,ID+.05,0);
            for clusternum = 1:size(outputcifti,2)
                clustercount_new = clustercount_new + 1;
                
                
                this_maxdice = 0;
                this_dicediff = 0;
                
                
                cluster_surfacearea = sum(surfacearea(logical(outputcifti(:,clusternum))));
                
                clusters_subs_IDs_dices_SAs_new(clustercount_new,:) = [s ID this_maxdice this_dicediff cluster_surfacearea];
                subnetworks_bycluster_new(logical(outputcifti(:,clusternum)),s) = clustercount_new;
            end
        end
    end
    disp(' ')
    
    save('System_clustering_new.mat','clusters_subs_IDs_dices_SAs_new','subnetworks_bycluster_new','clustercount','-v7.3')
else
    load('System_clustering_new.mat')
end

%%
load('System_clustering_old.mat')



%%
load('System_clustering_old.mat')

noldsubs = max(clusters_subs_IDs_dices_SAs(:,1));
clusters_subs_IDs_dices_SAs_new(:,1) = clusters_subs_IDs_dices_SAs_new(:,1) + noldsubs;
clusters_subs_IDs_dices_SAs = [clusters_subs_IDs_dices_SAs ; clusters_subs_IDs_dices_SAs_new];
subnetworks_bycluster_new(subnetworks_bycluster_new>0) = subnetworks_bycluster_new(subnetworks_bycluster_new>0) + clustercount;
subnetworks_bycluster = [subnetworks_bycluster subnetworks_bycluster_new];
clustercount = clustercount + clustercount_new;
nsubs = max(clusters_subs_IDs_dices_SAs(:,1));

distances = smartload('/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_standardsubcort_uint8.mat');
%distances = uint8(distances);

indices = {[]};
count = 0;

prevstring = [];
for s1 = 1:nsubs
    s1clus_dmat = [];
    for s2 = (max([(s1+1) (noldsubs+1)])):nsubs
        
       
            string = ['Matching subjects pairwise: Subject ' num2str(s1) ' vs Subject ' num2str(s2) ', '];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = [string 'Network #' num2str(max(IDs)) '  '];
            
            [s1matches,s2matches,s1clus_dmat] = networkClusters_match_preclustered(subnetworks_bycluster(:,s1),[find(clusters_subs_IDs_dices_SAs(:,1)==s1) clusters_subs_IDs_dices_SAs(clusters_subs_IDs_dices_SAs(:,1)==s1,2)],distances,s1clus_dmat,subnetworks_bycluster(:,s2),[find(clusters_subs_IDs_dices_SAs(:,1)==s2) clusters_subs_IDs_dices_SAs(clusters_subs_IDs_dices_SAs(:,1)==s2,2)],distthresh);
            
            for matchnum = 1:size(s1matches,2)
                
                s1_clusternums = unique(subnetworks_bycluster(logical(s1matches(:,matchnum)),s1));
                s2_clusternums = unique(subnetworks_bycluster(logical(s2matches(:,matchnum)),s2));
                
                count = count+1;
                indices{count,1} = s1_clusternums;
                indices{count,2} = s2_clusternums;
                
                %match_matrix(s1_clusternums,s2_clusternums) = 1;  match_matrix(s2_clusternums,s1_clusternums) = 1;
                
            end
        
    end
end
disp(' ')
clear distances

match_matrix_new = false(clustercount);
match_matrix_new(1:size(match_matrix,1),1:size(match_matrix,2)) = match_matrix;
clear match_matrix
for i = 1:size(indices,1)
    match_matrix_new(indices{i,1},indices{i,2}) = true;  match_matrix_new(indices{i,2},indices{i,1}) = true;
end

save('match_matrix_new.mat','match_matrix_new','-v7.3')

%%

mat2pajek_byindex(match_matrix_new,find(triu(match_matrix,1)),'pajek.net')
rawclrs = run_infomap_on_pajekfile('pajek.net',100);
dlmwrite(['rawassn.txt'],rawclrs,'\t')
delete('pajek.clu')
simple_assigns = modify_clrfile('simplify','rawassn.txt',2);
