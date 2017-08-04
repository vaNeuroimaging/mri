

MSCnames = {'MSC01','MSC02','MSC03','MSC04'};

ncortverts = 59412;
 
group_networks = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii'); 
networkIDs = unique(group_networks); networkIDs(networkIDs<1) = [];
groupclusters = zeros(size(group_networks,1),0);
groupclusternetworkIDs = [];
for ID = networkIDs(:)'
    clusters = metric_cluster_cifti(group_networks,ID-.5,ID+.5,20);
    groupclusters(:,end+1:end+size(clusters,2)) = clusters;
    groupclusternetworkIDs(end+1:end+size(clusters,2)) = ID;
end
groupclusters = groupclusters(1:ncortverts,:);
groupnetworks_noassign = find(group_networks < 1);



task_maps = {[1 2 5 7 9 11 13 15 17],[1],[1:2:23]};



sub_networks = cifti_read('Template_match_systems_Poldrome_MSC.dtseries.nii'); sub_networks = sub_networks(1:ncortverts,:);
sub_networks(groupnetworks_noassign,:) = 0;

subnet_similarity_mat = zeros(size(sub_networks,2));
subnet_similarity_mat_within = zeros(size(sub_networks,2));
subnet_similarity_mat_between = zeros(size(sub_networks,2));
groupnet_similarity_mat = zeros(size(sub_networks,2));
groupnet_similarity_mat_within = zeros(size(sub_networks,2));
groupnet_similarity_mat_between = zeros(size(sub_networks,2));

subact_similarity_mat = zeros(size(sub_networks,2));
groupact_similarity_mat = zeros(size(sub_networks,2));

for s1 = 1:size(sub_networks,2)
    
    if s1==1
        sub1timeseries = cifti_read('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_LR_timeseries.dtseries.nii');
        tmask = load('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_total_tmask.txt');
        sub1timeseries = sub1timeseries(1:ncortverts,logical(tmask));
    elseif s1 < size(sub_networks,2)
        MSCname = MSCnames{s1-1};
        tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
        [subjects tmasks] = textread(tmaskfile,'%s %s');
        for i = 1:length(subjects)
            if i == 1
                sub1timeseries = cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                tmask = load(tmasks{i});
            else
                sub1timeseries = [sub1timeseries cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'])];
                tmask = [tmask; load(tmasks{i})];
            end
        end
        sub1timeseries = sub1timeseries(1:ncortverts,logical(tmask));
        
        
        
    end
    
    groupcluster_timecourses = zeros(size(sub1timeseries,2),size(groupclusters,2));
    for cluster = 1:size(groupclusters,2)
        groupcluster_timecourses(:,cluster) = mean(sub1timeseries(logical(groupclusters(:,cluster)),:),1);
    end
    
    gclustermat_s1 = FisherTransform(paircorr_mod(groupcluster_timecourses));
    
    for s2 = 1:size(sub_networks,2)
        if s2 > s1
            
            MSCname = MSCnames{s2-1};
            tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
            [subjects tmasks] = textread(tmaskfile,'%s %s');
            for i = 1:length(subjects)
                if i == 1
                    sub2timeseries = cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                    tmask = load(tmasks{i});
                else
                    sub2timeseries = [sub2timeseries cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'])];
                    tmask = [tmask; load(tmasks{i})];
                end
            end
            sub2timeseries = sub2timeseries(1:ncortverts,logical(tmask));
            
            groupcluster_timecourses = zeros(size(sub2timeseries,2),size(groupclusters,2));
            for cluster = 1:size(groupclusters,2)
                groupcluster_timecourses(:,cluster) = mean(sub2timeseries(logical(groupclusters(:,cluster)),:),1);
            end
            
            gclustermat_s2 = FisherTransform(paircorr_mod(groupcluster_timecourses));
            
            triumat = logical(triu(ones(size(gclustermat_s1)),1));
            
            group_within_network_mat = logical(~logical(squareform(pdist(groupclusternetworkIDs'))) - eye(length(groupclusternetworkIDs)));
            group_between_network_mat = logical(squareform(pdist(groupclusternetworkIDs')));
            
            
            groupnet_similarity_mat(s1,s2) = FisherTransform(paircorr_mod(gclustermat_s1(triumat),gclustermat_s2(triumat)));
            groupnet_similarity_mat_within(s1,s2) = FisherTransform(paircorr_mod(gclustermat_s1(logical(group_within_network_mat.*triumat)),gclustermat_s2(logical(group_within_network_mat.*triumat))));
            groupnet_similarity_mat_between(s1,s2) = FisherTransform(paircorr_mod(gclustermat_s1(logical(group_between_network_mat.*triumat)),gclustermat_s2(logical(group_between_network_mat.*triumat))));
            
            
            disp(['Matching network clusters: subject ' num2str(s1) ' vs subject ' num2str(s2)'])
            [s1matches, s2matches] = networkClusters_match(sub_networks(:,s1),distances,sub_networks(:,s2));
            
            s1cluster_timecourses = zeros(size(sub1timeseries,2),size(s1matches,2));
            s2cluster_timecourses = zeros(size(sub2timeseries,2),size(s1matches,2));
            subclusternetworkIDs = zeros(size(s1matches,2));
            for cluster = 1:size(s1matches,2)
                s1cluster_timecourses(:,cluster) = mean(sub1timeseries(logical(s1matches(:,cluster)),:),1);
                s2cluster_timecourses(:,cluster) = mean(sub2timeseries(logical(s2matches(:,cluster)),:),1);
                subclusternetworkIDs(cluster,1) = mode(sub_networks(logical(s1matches(:,cluster)),s1));
            end
            
            sclustermat_s1 = FisherTransform(paircorr_mod(s1cluster_timecourses));
            sclustermat_s2 = FisherTransform(paircorr_mod(s2cluster_timecourses));
            triumat = logical(triu(ones(size(sclustermat_s1)),1));
            
            sub_within_network_mat = logical(~logical(squareform(pdist(subclusternetworkIDs))) - eye(length(subclusternetworkIDs)));
            sub_between_network_mat = logical(squareform(pdist(subclusternetworkIDs)));
            
            
            subnet_similarity_mat(s1,s2) = FisherTransform(paircorr_mod(sclustermat_s1(triumat),sclustermat_s2(triumat)));
            subnet_similarity_mat_within(s1,s2) = FisherTransform(paircorr_mod(sclustermat_s1(logical(sub_within_network_mat.*triumat)),sclustermat_s2(logical(sub_within_network_mat.*triumat))));
            subnet_similarity_mat_between(s1,s2) = FisherTransform(paircorr_mod(sclustermat_s1(logical(sub_between_network_mat.*triumat)),sclustermat_s2(logical(sub_between_network_mat.*triumat))));
            
            
            if s1 > 1
                
                groupcluster_activations_s1 = zeros(size(groupclusters,2),0);
                groupcluster_activations_s2 = zeros(size(groupclusters,2),0);
                subcluster_activations_s1 = zeros(size(s1matches,2),0);
                subcluster_activations_s2 = zeros(size(s1matches,2),0);
            
            MSCname = MSCnames{s1-1};
            task_files = {['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.motor_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.memory_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.mixed_taskmaps.dscalar.nii']};
            contrastcounter = 0;
            for filenum = 1:length(task_files)
                taskdata = cifti_read(task_files{filenum});
                for contrastnum = 1:length(task_maps{filenum})
                    contrast = task_maps{filenum}(contrastnum);
                    contrastcounter = contrastcounter +1;
                    
                    for cluster = 1:size(groupclusters,2)
                        groupcluster_activations_s1(cluster,contrastcounter) = mean(taskdata(logical(groupclusters(:,cluster)),contrast));
                    end
                    for cluster = 1:size(s1matches,2)
                        subcluster_activations_s1(cluster,contrastcounter) = mean(taskdata(logical(s1matches(:,cluster)),contrast));
                    end
                end
            end
            
            MSCname = MSCnames{s2-1};
            task_files = {['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.motor_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.memory_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.mixed_taskmaps.dscalar.nii']};
            contrastcounter = 0;
            for filenum = 1:length(task_files)
                taskdata = cifti_read(task_files{filenum});
                for contrastnum = 1:length(task_maps{filenum})
                    contrast = task_maps{filenum}(contrastnum);
                    contrastcounter = contrastcounter +1;
                    
                    for cluster = 1:size(groupclusters,2)
                        groupcluster_activations_s2(cluster,contrastcounter) = mean(taskdata(logical(groupclusters(:,cluster)),contrast));
                    end
                    for cluster = 1:size(s1matches,2)
                        subcluster_activations_s2(cluster,contrastcounter) = mean(taskdata(logical(s2matches(:,cluster)),contrast));
                    end
                end
            end
            
            groupact_similarity_mat(s1,s2) = FisherTransform(paircorr_mod(groupcluster_activations_s1(:),groupcluster_activations_s2(:)));
            subact_similarity_mat(s1,s2) = FisherTransform(paircorr_mod(subcluster_activations_s1(:),subcluster_activations_s2(:)));
            
            end
            
            
            
        end
    end
end

triumat = logical(triu(ones(size(subnet_similarity_mat)),1));

disp(['Mean group network cluster matrix REST similarities: ' num2str(mean(groupnet_similarity_mat(triumat))) ' +/- ' num2str(std(groupnet_similarity_mat(triumat)))])
disp(['Mean subject network cluster matrix REST similarities: ' num2str(mean(subnet_similarity_mat(triumat))) ' +/- ' num2str(std(subnet_similarity_mat(triumat)))])
[H,P,CI,STATS] = ttest(subnet_similarity_mat(triumat),groupnet_similarity_mat(triumat));
disp(['Paired t-test REST results: T=' num2str(STATS.tstat) ', p=' num2str(P)])

figure;imagesc(groupnet_similarity_mat,[.75 1]); colormap hot; title('Subject REST similarities: group network clusters'); colorbar
figure;imagesc(subnet_similarity_mat,[.75 1]); colormap hot; title('Subject REST similarities: subject network clusters'); colorbar

figure;imagesc(subnet_similarity_mat - groupnet_similarity_mat,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject REST similarities: subject network clusters > group network clusters'); colorbar

disp(['Mean group network cluster matrix REST similarities within network: ' num2str(mean(groupnet_similarity_mat_within(triumat))) ' +/- ' num2str(std(groupnet_similarity_mat_within(triumat)))])
disp(['Mean subject network cluster matrix REST similarities within network: ' num2str(mean(subnet_similarity_mat_within(triumat))) ' +/- ' num2str(std(subnet_similarity_mat_within(triumat)))])
[H,P,CI,STATS] = ttest(subnet_similarity_mat_within(triumat),groupnet_similarity_mat_within(triumat));
disp(['Paired t-test REST results: T=' num2str(STATS.tstat) ', p=' num2str(P)])

figure;imagesc(groupnet_similarity_mat_within,[.75 1]); colormap hot; title('Subject REST similarities within network: group network clusters'); colorbar
figure;imagesc(subnet_similarity_mat_within,[.75 1]); colormap hot; title('Subject REST similarities within network: subject network clusters'); colorbar

figure;imagesc(subnet_similarity_mat_within - groupnet_similarity_mat_within,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject REST similarities within network: subject network clusters > group network clusters'); colorbar

disp(['Mean group network cluster matrix REST similarities between network: ' num2str(mean(groupnet_similarity_mat_between(triumat))) ' +/- ' num2str(std(groupnet_similarity_mat_between(triumat)))])
disp(['Mean subject network cluster matrix REST similarities between network: ' num2str(mean(subnet_similarity_mat_between(triumat))) ' +/- ' num2str(std(subnet_similarity_mat_between(triumat)))])
[H,P,CI,STATS] = ttest(subnet_similarity_mat_between(triumat),groupnet_similarity_mat_between(triumat));
disp(['Paired t-test REST results: T=' num2str(STATS.tstat) ', p=' num2str(P)])

figure;imagesc(groupnet_similarity_mat_between,[.75 1]); colormap hot; title('Subject REST similarities between network: group network clusters'); colorbar
figure;imagesc(subnet_similarity_mat_between,[.75 1]); colormap hot; title('Subject REST similarities between network: subject network clusters'); colorbar

figure;imagesc(subnet_similarity_mat_between - groupnet_similarity_mat_between,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject REST similarities between network: subject network clusters > group network clusters'); colorbar

triumat(1,:) = 0;

disp(['Mean group network cluster matrix TASK similarities: ' num2str(mean(groupact_similarity_mat(triumat))) ' +/- ' num2str(std(groupact_similarity_mat(triumat)))])
disp(['Mean subject network cluster matrix TASK similarities: ' num2str(mean(subact_similarity_mat(triumat))) ' +/- ' num2str(std(subact_similarity_mat(triumat)))])
[H,P,CI,STATS] = ttest(subact_similarity_mat(triumat),groupact_similarity_mat(triumat));
disp(['Paired t-test TASK results: T=' num2str(STATS.tstat) ', p=' num2str(P)])

figure;imagesc(groupact_similarity_mat,[.75 1]); colormap hot; title('Subject TASK similarities: group network clusters'); colorbar
figure;imagesc(subact_similarity_mat,[.75 1]); colormap hot; title('Subject TASK similarities: subject network clusters'); colorbar

figure;imagesc(subact_similarity_mat - groupact_similarity_mat,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject TASK similarities: subject network clusters > group network clusters'); colorbar
            
            