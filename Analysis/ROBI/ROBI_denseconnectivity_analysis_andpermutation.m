subjects = {'Arsenio'};%{'ROBI002','ROBI004','ROBI005',
outfolder = '/home/data/Analysis/ROBI';
postsessions = [101 102 103];
bootstrap_iterations = 100;

height_threshs = [.08 : .01 : .15];
alpha = .05;

ncortverts = 59412;

write_dconns = 1;
do_bootstrap = 0;

%%

for s = 1:length(subjects)
    
    %% load data
    
    datafile = ['/home/data/subjects/' subjects{s} '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = ft_read_cifti_mod(datafile);
    
    datastruct = data; datastruct.data = [];
    
    runs_sessionsfile = ['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_runs_sessions.txt'];
    tmask = load(['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_all_tmask.txt']);
    
    runs_sessions = load(runs_sessionsfile);
    runs_sessions = runs_sessions(logical(tmask),:);
    sessions = runs_sessions(:,2);
    runs = runs_sessions(:,1);
    clear run_sessions
    
    
    
    
    %% pre vs post
    for postsess = postsessions
        
        corr_pre = single(paircorr_mod(data.data(:,sessions<100)'));
        corr_pre(isnan(corr_pre)) = 0;
        corr_pre = FisherTransform(corr_pre);
        
        
        
        
        if any(sessions==postsess)
            postname = ['post' num2str(postsess)];
            
            corr_post = single(paircorr_mod(data.data(:,sessions==postsess)'));
            corr_post(isnan(corr_post)) = 0;
            corr_post = FisherTransform(corr_post);
            
            
            corr_postminpre = corr_post - corr_pre;
            
            if write_dconns
                datastruct.data = corr_pre;
                clear corr_pre
                datastruct.dimord = 'pos_pos';
                ft_write_cifti_mod([outfolder '/' subjects{s} '_pre'],datastruct);
                datastruct.data = [];
                
                
                datastruct.data = corr_post;
                clear corr_post
                ft_write_cifti_mod([outfolder '/' subjects{s} '_' postname],datastruct);
                datastruct.data = [];
                
                datastruct.data = corr_postminpre;
                ft_write_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre'],datastruct);
                datastruct.data = [];
            end
            
            clear corr_pre corr_post
            
            postminpre_mean = mean(abs(corr_postminpre),2);
            clear corr_postminpre
            datastruct.data = postminpre_mean;
            datastruct.dimord = 'pos_time';
            ft_write_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre'],datastruct);
            
        end
    end
    
    
    
    if(do_bootstrap)
    %% bootstrap
    for postsess = postsessions
        if any(sessions==postsess)
            postname = ['post' num2str(postsess)];
            
            
            rand_postminpre = zeros(size(data.data,1),bootstrap_iterations);
            
            prerunIDs = unique(runs(sessions<100));
            postrunIDs = unique(runs(sessions==postsess));
            runIDs = [prerunIDs; postrunIDs];
            runIDcomparisons = [ones(size(prerunIDs)); (ones(size(postrunIDs)).*2)];
            
            for iter = 1:bootstrap_iterations
                disp([subjects{s} ': bootstrap permutation ' num2str(iter)])
                
                randcomparison = runIDcomparisons(randperm(length(runIDcomparisons)));
                randprepostindex = zeros(size(sessions));
                for runnum = 1:length(runIDs)
                    randprepostindex(runs==runIDs(runnum)) = randcomparison(runnum);
                end
                
                corr_pre = single(paircorr_mod(data.data(:,randprepostindex==1)'));
                if iter==1
                    nanmat = isnan(corr_pre);
                end
                corr_pre(nanmat) = 0;
                corr_pre = FisherTransform(corr_pre);
                
                corr_post = single(paircorr_mod(data.data(:,randprepostindex==2)'));
                corr_post(nanmat) = 0;
                corr_post = FisherTransform(corr_post);
                
                rand_postminpre(:,iter) = mean(abs(corr_post - corr_pre),2);
                
                clear corr_pre corr_post
            end
            datastruct.data = rand_postminpre;
            datastruct.dimord = 'pos_time';
            ft_write_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre_random'],datastruct);
            
             
            %% z-stat
            
            datastruct = ft_read_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre.dtseries.nii']);
            postminpre_mean = datastruct.data;
            
            zstat = (postminpre_mean - mean(rand_postminpre,2)) ./ std(rand_postminpre,[],2);
            
            datastruct.data = zstat;
            ft_write_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre_zstatVrandom'],datastruct);
            
            
            %% cluster correct
            
            default_surfaceareas = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_cortexonly.dtseries.nii';
            surfaceareas = ft_read_cifti_mod(default_surfaceareas);
            
            neighbors = cifti_neighbors([outfolder '/' subjects{s} '_' postname 'minpre.dtseries.nii']);
            
            datastruct = ft_read_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre_random.dtseries.nii']);
            rand_postminpre = datastruct.data;
            
            datastruct = ft_read_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre.dtseries.nii']);
            postminpre_mean = datastruct.data;
            
            postminpre_mean_cluscorrected = zeros(size(postminpre_mean,1),length(height_threshs));
            datastruct.data = [];
            datastruct.dimord = 'scalar_pos';
            datastruct.mapname = cell(1,length(height_threshs));
            
            for height_thresh_num = 1:length(height_threshs)
                height_thresh = height_threshs(height_thresh_num);
                max_random_clustersizes = zeros(bootstrap_iterations,1);
                for iter = 1:bootstrap_iterations
                    disp(['Clustering iteration ' num2str(iter)])
                    permuted_clustereddeltas = cifti_cluster_surface(rand_postminpre(1:ncortverts,iter),height_thresh,inf,0,surfaceareas.data,neighbors);
                    permuted_clustersizes = zeros(size(permuted_clustereddeltas,2),1);
                    for i = 1:size(permuted_clustereddeltas,2)
                        permuted_clustersizes(i) = sum(permuted_clustereddeltas(:,i) .* surfaceareas.data,1);
                    end
                    max_random_clustersizes(iter) = max(permuted_clustersizes);
                    sorted_max_random_clustersizes = sort(max_random_clustersizes,'descend');
                end
                
                clustersize_thresh = sorted_max_random_clustersizes(round(bootstrap_iterations * alpha));
                disp(['Height threshold = ' num2str(height_thresh) ', Corrected alpha = ' num2str(alpha) ': cluster size threshold = ' num2str(clustersize_thresh) 'mm'])
                
                clustereddeltas = cifti_cluster_surface(postminpre_mean(1:ncortverts),height_thresh,inf,clustersize_thresh,surfaceareas.data,neighbors);
                clustereddeltas(ncortverts+1:size(postminpre_mean,1),:) = 0;
                output = postminpre_mean .* sum(clustereddeltas,2);
                
                postminpre_mean_cluscorrected(:,height_thresh_num) = output;
                datastruct.mapname{1,height_thresh_num} = ['Height = ' num2str(height_thresh) ', Alpha = ' num2str(alpha) ': k = ' num2str(clustersize_thresh) 'mm'];
            end
            datastruct.data = postminpre_mean_cluscorrected;
            ft_write_cifti_mod([outfolder '/' subjects{s} '_' postname 'minpre_clustercorrected'],datastruct);
            
        end
    end
    end
    
    
    
end