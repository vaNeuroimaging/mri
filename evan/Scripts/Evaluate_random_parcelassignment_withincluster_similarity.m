

MSCnames = {'MSC01','MSC02','MSC03','MSC04'};

ncortverts = 59412;

xdistance = 15;

nrandomiter = 100;

%distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
distances = distances(1:ncortverts,1:ncortverts);

subparcelfile = ['Parcels_Poldrome_MSC.dtseries.nii'];
subparcels = ft_read_cifti_mod(subparcelfile); subparcels = subparcels.data(1:ncortverts,:);
sub_networks = ft_read_cifti_mod('Template_match_systems_Poldrome_MSC.dtseries.nii'); sub_networks = sub_networks.data(1:ncortverts,:);

assignmentsfilename = 'Assignments_xd15_1.0thr_merged_adjacent.dtseries.nii';
%assigned_parcels = cifti_read(assignmentsfilename);
assigned_parcels = ft_read_cifti_mod(assignmentsfilename); assigned_parcels = assigned_parcels.data;
assigned_parcelIDs = unique(assigned_parcels); assigned_parcelIDs(assigned_parcelIDs<1) = [];

assigned_parcel_networks = zeros(length(assigned_parcelIDs),1);
for i = 1:length(assigned_parcelIDs)
    IDs_byvert = sub_networks(assigned_parcels(1:ncortverts,:)==assigned_parcelIDs(i));
    IDs_byvert(IDs_byvert==0) = [];
    assigned_parcel_networks(i) = mode(IDs_byvert);
end

%group_parcels = cifti_read('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii');
group_parcels = ft_read_cifti_mod('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii'); group_parcels = group_parcels.data;
group_parcelIDs = unique(group_parcels); group_parcelIDs(group_parcelIDs==0) = [];

task_maps = {[1 2 5 7 9 11 13 15 17],[1],[1:2:23]};

parcels_subs_IDs = zeros(0,2);
for s = 1:size(subparcels,2);
    IDs = unique(subparcels(:,s));
    IDs(IDs==0) = [];
    parcels_subs_IDs = [parcels_subs_IDs ; [repmat(s,length(IDs),1) IDs]];
end
parcels_subs_IDs(:,end+1) = 0;

parcel_centroids = zeros(size(parcels_subs_IDs,1),1);

for i = 1:size(parcels_subs_IDs,1)
    thissub = parcels_subs_IDs(i,1);
    thisID = parcels_subs_IDs(i,2);
    parcelinds = find(subparcels(:,thissub)==thisID);
    networkIDs_byvert = sub_networks(parcelinds,thissub);
    networkIDs_byvert(networkIDs_byvert==0) = [];
    thisnetworkID = mode(networkIDs_byvert);
    parcels_subs_IDs(i,3) = thisnetworkID;
    
    parcel_totaldistances = sum(distances(parcelinds,parcelinds),2);
    [ign mini] = min(parcel_totaldistances);
    parcel_centroids(i) = parcelinds(mini);
    
end


group_restsimilaritymat = zeros(size(sub_networks,2));
assigned_restsimilaritymat = zeros(size(sub_networks,2));
closestdistance_assigned_restsimilaritymat = zeros(size(sub_networks,2));
withincluster_assigned_restsimilaritymat = zeros(size(sub_networks,2),size(sub_networks,2),nrandomiter);
withincluster_indistance_assigned_restsimilaritymat = zeros(size(sub_networks,2),size(sub_networks,2),nrandomiter);

group_tasksimilaritymat = zeros(size(sub_networks,2));
assigned_tasksimilaritymat = zeros(size(sub_networks,2));
closestdistance_assigned_tasksimilaritymat = zeros(size(sub_networks,2));
withincluster_assigned_tasksimilaritymat = zeros(size(sub_networks,2),size(sub_networks,2),nrandomiter);
withincluster_indistance_assigned_tasksimilaritymat = zeros(size(sub_networks,2),size(sub_networks,2),nrandomiter);

% assigned_restsimilarity_mat = zeros(size(sub_networks,2));
% assigned_restsimilarity_mat_within = zeros(size(sub_networks,2));
% assigned_restsimilarity_mat_between = zeros(size(sub_networks,2));
% groupnet_similarity_mat = zeros(size(sub_networks,2));
% groupnet_similarity_mat_within = zeros(size(sub_networks,2));
% groupnet_similarity_mat_between = zeros(size(sub_networks,2));
% 
% subact_similarity_mat = zeros(size(sub_networks,2));
% groupact_similarity_mat = zeros(size(sub_networks,2));
for i = 1:size(sub_networks,2)
    subparcel_correlmat{i} = [];
end

for s1 = 1:size(sub_networks,2)
    
    if s1==1
        %sub1timeseries = cifti_read('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_LR_timeseries.dtseries.nii');
        sub1timeseries = ft_read_cifti_mod('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_LR_timeseries.dtseries.nii'); sub1timeseries = sub1timeseries.data;
        tmask = load('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_total_tmask.txt');
        sub1timeseries = sub1timeseries(1:ncortverts,logical(tmask));
    elseif s1 < size(sub_networks,2)
        MSCname = MSCnames{s1-1};
        tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
        [subjects tmasks] = textread(tmaskfile,'%s %s');
        for i = 1:length(subjects)
            if i == 1
                sub1timeseries = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                sub1timeseries = sub1timeseries.data;
                tmask = load(tmasks{i});
            else
                temp = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                sub1timeseries = [sub1timeseries temp.data];
                clear temp
                tmask = [tmask; load(tmasks{i})];
            end
        end
        sub1timeseries = sub1timeseries(1:ncortverts,logical(tmask));
    end
    
    sub1parcelinds = find(parcels_subs_IDs(:,1)==s1);
    sub1parcelIDs = parcels_subs_IDs(sub1parcelinds,2);
    
    if isempty(subparcel_correlmat{s1})
        timecourses = zeros(size(sub1timeseries,2),length(sub1parcelIDs));
        for i = 1:length(sub1parcelIDs)
            timecourses(:,i) = mean(sub1timeseries(subparcels(:,s1)==sub1parcelIDs(i),:),1);
        end
        subparcel_correlmat{s1} = FisherTransform(paircorr_mod(timecourses));
        
        
        timecourses = zeros(size(sub1timeseries,2),length(assigned_parcelIDs));
        for i = 1:length(assigned_parcelIDs)
            if any(assigned_parcels(:,s1)==assigned_parcelIDs(i))
                timecourses(:,i) = mean(sub1timeseries(assigned_parcels(:,s1)==assigned_parcelIDs(i),:),1);
                
            end
        end
        assignedparcel_correlmat{s1} = FisherTransform(paircorr_mod(timecourses));
        
        timecourses = zeros(size(sub1timeseries,2),length(group_parcelIDs));
        for i = 1:length(group_parcelIDs)
            if any(group_parcels==group_parcelIDs(i))
                timecourses(:,i) = mean(sub1timeseries(group_parcels==group_parcelIDs(i),:),1);
            end
        end
        groupparcel_correlmat{s1} = FisherTransform(paircorr_mod(timecourses));
        
    end
        
        
    

    
    for s2 = 1:size(sub_networks,2)
        if s2 > s1
            
            MSCname = MSCnames{s2-1};
            tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
            [subjects tmasks] = textread(tmaskfile,'%s %s');
            for i = 1:length(subjects)
                if i == 1
                    sub2timeseries = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                    sub2timeseries = sub2timeseries.data;
                    tmask = load(tmasks{i});
                else
                    temp = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                    sub2timeseries = [sub2timeseries temp.data]; 
                    clear temp
                    tmask = [tmask; load(tmasks{i})];
                end
            end
            sub2timeseries = sub2timeseries(1:ncortverts,logical(tmask));
            
            sub2parcelinds = find(parcels_subs_IDs(:,1)==s2);
            sub2parcelIDs = parcels_subs_IDs(sub2parcelinds,2);
            
            if isempty(subparcel_correlmat{s2})
                timecourses = zeros(size(sub2timeseries,2),length(sub2parcelIDs));
                for i = 1:length(sub2parcelIDs)
                    timecourses(:,i) = mean(sub2timeseries(subparcels(:,s2)==sub2parcelIDs(i),:),1);
                end
                subparcel_correlmat{s2} = FisherTransform(paircorr_mod(timecourses));
                
                timecourses = zeros(size(sub2timeseries,2),length(assigned_parcelIDs));
                for i = 1:length(assigned_parcelIDs)
                    if any(assigned_parcels(:,s2)==assigned_parcelIDs(i))
                        timecourses(:,i) = mean(sub2timeseries(assigned_parcels(:,s2)==assigned_parcelIDs(i),:),1);
                    end
                end
                assignedparcel_correlmat{s2} = FisherTransform(paircorr_mod(timecourses));
                
                timecourses = zeros(size(sub2timeseries,2),length(group_parcelIDs));
                for i = 1:length(group_parcelIDs)
                    if any(group_parcels==group_parcelIDs(i))
                        timecourses(:,i) = mean(sub2timeseries(group_parcels==group_parcelIDs(i),:),1);
                    end
                end
                groupparcel_correlmat{s2} = FisherTransform(paircorr_mod(timecourses));
            end
            
            triumat = triu(true(size(groupparcel_correlmat{s1})),1);
            group_restsimilaritymat(s1,s2) = FisherTransform(paircorr_mod(groupparcel_correlmat{s1}(triumat),groupparcel_correlmat{s2}(triumat)));
            
            
            
            s1_parcelIDs = unique(assigned_parcels(:,s1)); s1_parcelIDs(s1_parcelIDs<1) = [];
            s2_parcelIDs = unique(assigned_parcels(:,s2)); s2_parcelIDs(s2_parcelIDs<1) = [];
            commonassignedIDs = intersect(s1_parcelIDs,s2_parcelIDs);
            [ign commonassignedID_inds ign2] = intersect(assigned_parcelIDs,commonassignedIDs);
            
            triumat = triu(true(length(commonassignedID_inds)),1);
            s1mat = assignedparcel_correlmat{s1}(commonassignedID_inds,commonassignedID_inds);
            s2mat = assignedparcel_correlmat{s2}(commonassignedID_inds,commonassignedID_inds);
            assigned_restsimilaritymat(s1,s2) = FisherTransform(paircorr_mod(s1mat(triumat),s2mat(triumat)));
            
            if ((s1==1) && (s2==2)) || ((s1==3) && (s2==4)) || ((s1==1) && (s2==5))
                1;
            end
            
            
            
            
            full_distmat = distances(parcel_centroids(sub1parcelinds),parcel_centroids(sub2parcelinds));
            full_distmat(full_distmat > xdistance) = Inf;
            [assign_bydistance, ign] = munkres(full_distmat);
            assign_bydistance_inds = find(assign_bydistance);
            
            s1_parcelmat = subparcel_correlmat{s1}(assign_bydistance_inds,assign_bydistance_inds);
            s2_parcelmat = subparcel_correlmat{s2}(assign_bydistance(assign_bydistance_inds),assign_bydistance(assign_bydistance_inds));
            triumat = logical(triu(ones(size(s1_parcelmat)),1));
            
            closestdistance_assigned_restsimilaritymat(s1,s2) = FisherTransform(paircorr_mod(s1_parcelmat(triumat),s2_parcelmat(triumat)));
            
            
            
            
            
            disp(['Matching network clusters: subject ' num2str(s1) ' vs subject ' num2str(s2)'])
            [s1matches, s2matches] = networkClusters_match(sub_networks(:,s1),distances,sub_networks(:,s2));
            
            
            sub1parcels_incluster = zeros(length(sub1parcelIDs),1);
            for i = 1:length(sub1parcelIDs)
            
                numparcelverts = nnz(subparcels(:,s1)==sub1parcelIDs(i));
                
                numparcelverts_incluster = zeros(size(s1matches,2),1);
                for clusternum = 1:size(s1matches,2)
                    numparcelverts_incluster(clusternum) = nnz((subparcels(:,s1)==sub1parcelIDs(i)) .* s1matches(:,clusternum));
                end
                [maxverts_inclus maxi] = max(numparcelverts_incluster);
                
                [modeval nummode] = mode(sub_networks(subparcels(:,s1)==sub1parcelIDs(i),s1));
                
                if maxverts_inclus == nummode
                    
                    sub1parcels_incluster(i) = maxi;
                    
                end
            end
            
            sub2parcels_incluster = zeros(length(sub2parcelIDs),1);
            for i = 1:length(sub2parcelIDs)
                
                numparcelverts = nnz(subparcels(:,s2)==sub2parcelIDs(i));
                
                numparcelverts_incluster = zeros(size(s2matches,2),1);
                for clusternum = 1:size(s2matches,2)
                    numparcelverts_incluster(clusternum) = nnz((subparcels(:,s2)==sub2parcelIDs(i)) .* s2matches(:,clusternum));
                end
                [maxverts_inclus maxi] = max(numparcelverts_incluster);
                
                [modeval nummode] = mode(sub_networks(subparcels(:,s2)==sub2parcelIDs(i),s2));
                
                if maxverts_inclus == nummode
                    
                    sub2parcels_incluster(i) = maxi;
                    
                end
            end
            
            assign_bycluster = zeros(length(sub1parcelIDs),nrandomiter);
            assign_bycluster_indistance = zeros(length(sub1parcelIDs),nrandomiter);
            
            for iter = 1:nrandomiter
                
                
%                 for clusternum = 1:size(s2matches,2)
%                     sub1parcelinds_inthiscluster = find(sub1parcels_incluster==clusternum);
%                     sub2parcelinds_inthiscluster = find(sub2parcels_incluster==clusternum);
%                     randorder_s1 = randperm(length(sub1parcelinds_inthiscluster));
%                     randorder_s2 = randperm(length(sub2parcelinds_inthiscluster));
%                     
%                     minlength = min(length(sub1parcelinds_inthiscluster),length(sub2parcelinds_inthiscluster));
%                     
%                     assign_bycluster(sub1parcelinds_inthiscluster(randorder_s1(1:minlength)),iter) = sub2parcelinds_inthiscluster(randorder_s2(1:minlength));
%                     
%                 end
                
                %differentclustermat = logical(pdist2(sub1parcels_incluster,sub2parcels_incluster));
                randmat = rand(length(sub1parcelIDs),length(sub2parcelIDs));
                randmat_distconstrained = rand(length(sub1parcelIDs),length(sub2parcelIDs));
                randmat_distconstrained(isinf(full_distmat)) = Inf;
                
                for clusternum = 1:size(s2matches,2)
                    sub1parcelinds_inthiscluster = find(sub1parcels_incluster==clusternum);
                    sub2parcelinds_inthiscluster = find(sub2parcels_incluster==clusternum);
                    
                    assignedinds = munkres(randmat(sub1parcelinds_inthiscluster,sub2parcelinds_inthiscluster));
                    assign_bycluster(sub1parcelinds_inthiscluster(assignedinds>0),iter) = sub2parcelinds_inthiscluster(assignedinds(assignedinds>0));
                    
                    assignedinds = munkres(randmat_distconstrained(sub1parcelinds_inthiscluster,sub2parcelinds_inthiscluster));
                    assign_bycluster_indistance(sub1parcelinds_inthiscluster(assignedinds>0),iter) = sub2parcelinds_inthiscluster(assignedinds(assignedinds>0));
                end
                
                assign_bycluster_inds = find(assign_bycluster(:,iter));
                
                s1_parcelmat = subparcel_correlmat{s1}(assign_bycluster_inds,assign_bycluster_inds);
                s2_parcelmat = subparcel_correlmat{s2}(assign_bycluster(assign_bycluster_inds,iter),assign_bycluster(assign_bycluster_inds,iter));
                triumat = logical(triu(ones(size(s1_parcelmat)),1));
                
                withincluster_assigned_restsimilaritymat(s1,s2,iter) = FisherTransform(paircorr_mod(s1_parcelmat(triumat),s2_parcelmat(triumat)));
                
                
                assign_bycluster_indistance_inds = find(assign_bycluster_indistance(:,iter));
                
                s1_parcelmat = subparcel_correlmat{s1}(assign_bycluster_indistance_inds,assign_bycluster_indistance_inds);
                s2_parcelmat = subparcel_correlmat{s2}(assign_bycluster_indistance(assign_bycluster_indistance_inds,iter),assign_bycluster_indistance(assign_bycluster_indistance_inds,iter));
                triumat = logical(triu(ones(size(s1_parcelmat)),1));
                
                withincluster_indistance_assigned_restsimilaritymat(s1,s2,iter) = FisherTransform(paircorr_mod(s1_parcelmat(triumat),s2_parcelmat(triumat)));
            end
               
                
            
            
            
            
            
            
            
            
            
            
            
            
            if s1 > 1
                
                %                 groupcluster_activations_s1 = zeros(size(groupclusters,2),0);
                %                 groupcluster_activations_s2 = zeros(size(groupclusters,2),0);
                %                 subcluster_activations_s1 = zeros(size(s1matches,2),0);
                %                 subcluster_activations_s2 = zeros(size(s1matches,2),0);
                
                subparcel_activations_s1 = zeros(length(sub1parcelIDs),0);
                subparcel_activations_s2 = zeros(length(sub2parcelIDs),0);
                
                assignedparcel_activations_s1 = zeros(length(commonassignedIDs),0);
                assignedparcel_activations_s2 = zeros(length(commonassignedIDs),0);
                
                groupparcel_activations_s1 = zeros(length(group_parcelIDs),0);
                groupparcel_activations_s2 = zeros(length(group_parcelIDs),0);
                
                
                MSCname = MSCnames{s1-1};
                task_files = {['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.motor_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.memory_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.mixed_taskmaps.dscalar.nii']};
                contrastcounter = 0;
                for filenum = 1:length(task_files)
                    taskdata = ft_read_cifti_mod(task_files{filenum}); taskdata = taskdata.data;
                    for contrastnum = 1:length(task_maps{filenum})
                        contrast = task_maps{filenum}(contrastnum);
                        contrastcounter = contrastcounter +1;
                        
                        for parcelnum = 1:length(sub1parcelIDs)
                            subparcel_activations_s1(parcelnum,contrastcounter) = mean(taskdata(subparcels(:,s1)==sub1parcelIDs(parcelnum),contrast));
                        end
                        for parcelnum = 1:length(commonassignedIDs)
                            assignedparcel_activations_s1(parcelnum,contrastcounter) = mean(taskdata(assigned_parcels(:,s1)==commonassignedIDs(parcelnum),contrast));
                        end
                        for parcelnum = 1:length(group_parcelIDs)
                            groupparcel_activations_s1(parcelnum,contrastcounter) = mean(taskdata(group_parcels==group_parcelIDs(parcelnum),contrast));
                        end
                        
                    end
                end
                
                MSCname = MSCnames{s2-1};
                task_files = {['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.motor_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.memory_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.mixed_taskmaps.dscalar.nii']};
                contrastcounter = 0;
                for filenum = 1:length(task_files)
                    taskdata = ft_read_cifti_mod(task_files{filenum}); taskdata = taskdata.data;
                    for contrastnum = 1:length(task_maps{filenum})
                        contrast = task_maps{filenum}(contrastnum);
                        contrastcounter = contrastcounter +1;
                        
                        for parcelnum = 1:length(sub2parcelIDs)
                            subparcel_activations_s2(parcelnum,contrastcounter) = mean(taskdata(subparcels(:,s2)==sub2parcelIDs(parcelnum),contrast));
                        end
                        for parcelnum = 1:length(commonassignedIDs)
                            assignedparcel_activations_s2(parcelnum,contrastcounter) = mean(taskdata(assigned_parcels(:,s2)==commonassignedIDs(parcelnum),contrast));
                        end
                        for parcelnum = 1:length(group_parcelIDs)
                            groupparcel_activations_s2(parcelnum,contrastcounter) = mean(taskdata(group_parcels==group_parcelIDs(parcelnum),contrast));
                        end
                    end
                end
                
                group_tasksimilaritymat(s1,s2) = FisherTransform(paircorr_mod(groupparcel_activations_s1(:),groupparcel_activations_s2(:)));
                assigned_tasksimilaritymat(s1,s2) = FisherTransform(paircorr_mod(assignedparcel_activations_s1(:),assignedparcel_activations_s2(:)));
                
                s1_activations_bydistance = subparcel_activations_s1(assign_bydistance_inds,:);
                s2_activations_bydistance = subparcel_activations_s2(assign_bydistance(assign_bydistance_inds),:);
                closestdistance_assigned_tasksimilaritymat(s1,s2) = FisherTransform(paircorr_mod(s1_activations_bydistance(:),s2_activations_bydistance(:)));
                
                for iter = 1:nrandomiter
                    assign_bycluster_inds = find(assign_bycluster(:,iter));
                    s1_activations_withincluster = subparcel_activations_s1(assign_bycluster_inds,:);
                    s2_activations_withincluster = subparcel_activations_s2(assign_bycluster(assign_bycluster_inds,iter),:);
                    
                    withincluster_assigned_tasksimilaritymat(s1,s2,iter) = FisherTransform(paircorr_mod(s1_activations_withincluster(:),s2_activations_withincluster(:)));
                    
                    
                    assign_bycluster_indistance_inds = find(assign_bycluster_indistance(:,iter));
                    s1_activations_withincluster = subparcel_activations_s1(assign_bycluster_indistance_inds,:);
                    s2_activations_withincluster = subparcel_activations_s2(assign_bycluster_indistance(assign_bycluster_indistance_inds,iter),:);
                    
                    withincluster_indistance_assigned_tasksimilaritymat(s1,s2,iter) = FisherTransform(paircorr_mod(s1_activations_withincluster(:),s2_activations_withincluster(:)));
                    
                end
            end
            
            
            
        end
    end
end

avg_withincluster_assigned_restsimilaritymat = mean(withincluster_assigned_restsimilaritymat,3);
avg_withincluster_assigned_tasksimilaritymat = mean(withincluster_assigned_tasksimilaritymat,3);

avg_withincluster_indistance_assigned_restsimilaritymat = mean(withincluster_indistance_assigned_restsimilaritymat,3);
avg_withincluster_indistance_assigned_tasksimilaritymat = mean(withincluster_indistance_assigned_tasksimilaritymat,3);

triumat = logical(triu(ones(size(assigned_restsimilaritymat)),1));

disp(['Mean group parcel REST similarities: ' num2str(mean(group_restsimilaritymat(triumat))) ' +/- ' num2str(std(group_restsimilaritymat(triumat)))])
disp(['Mean assigned matrix REST similarities: ' num2str(mean(assigned_restsimilaritymat(triumat))) ' +/- ' num2str(std(assigned_restsimilaritymat(triumat)))])
disp(['Mean mindist parcel REST similarities: ' num2str(mean(closestdistance_assigned_restsimilaritymat(triumat))) ' +/- ' num2str(std(closestdistance_assigned_restsimilaritymat(triumat)))])
disp(['Mean rand within cluster matrix REST similarities: ' num2str(mean(avg_withincluster_assigned_restsimilaritymat(triumat))) ' +/- ' num2str(std(avg_withincluster_assigned_restsimilaritymat(triumat)))])
disp(['Mean rand within cluster within distance matrix REST similarities: ' num2str(mean(avg_withincluster_indistance_assigned_restsimilaritymat(triumat))) ' +/- ' num2str(std(avg_withincluster_indistance_assigned_restsimilaritymat(triumat)))])


figure;imagesc(group_restsimilaritymat,[.5 .9]); colormap hot; title('Subject REST similarities: group parcels'); colorbar
figure;imagesc(assigned_restsimilaritymat,[.5 .9]); colormap hot; title('Subject REST similarities: assigned parcels'); colorbar

figure;imagesc(assigned_restsimilaritymat - group_restsimilaritymat,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject REST similarities: assigned parcels > group parcels'); colorbar



triumat_task = triumat; triumat_task(1,:) = 0;
disp(['Mean group parcel TASK similarities: ' num2str(mean(group_tasksimilaritymat(triumat_task))) ' +/- ' num2str(std(group_tasksimilaritymat(triumat_task)))])
disp(['Mean assigned matrix TASK similarities: ' num2str(mean(assigned_tasksimilaritymat(triumat_task))) ' +/- ' num2str(std(assigned_tasksimilaritymat(triumat_task)))])
disp(['Mean mindist parcel TASK similarities: ' num2str(mean(closestdistance_assigned_tasksimilaritymat(triumat_task))) ' +/- ' num2str(std(closestdistance_assigned_tasksimilaritymat(triumat_task)))])
disp(['Mean rand within cluster matrix TASK similarities: ' num2str(mean(avg_withincluster_assigned_tasksimilaritymat(triumat_task))) ' +/- ' num2str(std(avg_withincluster_assigned_tasksimilaritymat(triumat_task)))])
disp(['Mean rand within cluster within distance matrix TASK similarities: ' num2str(mean(avg_withincluster_indistance_assigned_tasksimilaritymat(triumat_task))) ' +/- ' num2str(std(avg_withincluster_indistance_assigned_tasksimilaritymat(triumat_task)))])


figure;imagesc(group_tasksimilaritymat,[.5 .9]); colormap hot; title('Subject TASK similarities: group parcels'); colorbar
figure;imagesc(assigned_tasksimilaritymat,[.5 .9]); colormap hot; title('Subject TASK similarities: assigned parcels'); colorbar

figure;imagesc(assigned_tasksimilaritymat - group_tasksimilaritymat,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject TASK similarities: assigned parcels > group parcels'); colorbar


%%
%[H,P,CI,STATS] = ttest(assigned_restsimilarity_mat(triumat),group_restsimilaritymat(triumat));
%disp(['Paired t-test REST results: T=' num2str(STATS.tstat) ', p=' num2str(P)])

%figure;imagesc(group_restsimilaritymat,[.75 1]); colormap hot; title('Subject REST similarities: group network clusters'); colorbar
%figure;imagesc(assigned_restsimilarity_mat,[.75 1]); colormap hot; title('Subject REST similarities: subject network clusters'); colorbar
% 
% figure;imagesc(assigned_restsimilarity_mat - group_restsimilaritymat,[-.2 .2]); 
% colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
% hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
% combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
% colormap(combined);
% title('Subject REST similarities: subject network clusters > group network clusters'); colorbar

outputfilename = 'Similarities.txt';
delete([outputfilename]);
fid = fopen([outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Assignment','State','SubjectComparison','Similarity'); %write the output file header
fclose(fid);
dlmwrite([outputfilename],' ','-append');

for s1 = 1:size(sub_networks,2)
    for s2 = 1:size(sub_networks,2)
        if s2 > s1
            
            texttowrite = ['Group   Rest   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(group_restsimilaritymat(s1,s2))];  %save the data as a string to be written to the output
            dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
            texttowrite = ['Assigned   Rest   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(assigned_restsimilaritymat(s1,s2))];  %save the data as a string to be written to the output
            dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
            texttowrite = ['Closest   Rest   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(closestdistance_assigned_restsimilaritymat(s1,s2))];  %save the data as a string to be written to the output
            dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
            texttowrite = ['Rand_withinCluster   Rest   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(avg_withincluster_assigned_restsimilaritymat(s1,s2))];  %save the data as a string to be written to the output
            dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
            texttowrite = ['Rand_withinCluster_close   Rest   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(avg_withincluster_indistance_assigned_restsimilaritymat(s1,s2))];  %save the data as a string to be written to the output
            dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
            
            if s1 > 1
                texttowrite = ['Group   Task   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(group_tasksimilaritymat(s1,s2))];  %save the data as a string to be written to the output
                dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
                texttowrite = ['Assigned   Task   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(assigned_tasksimilaritymat(s1,s2))];  %save the data as a string to be written to the output
                dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
                texttowrite = ['Closest   Task   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(closestdistance_assigned_tasksimilaritymat(s1,s2))];  %save the data as a string to be written to the output
                dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
                texttowrite = ['Rand_withinCluster   Task   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(avg_withincluster_assigned_tasksimilaritymat(s1,s2))];  %save the data as a string to be written to the output
                dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
                texttowrite = ['Rand_withinCluster_close   Task   Sub' num2str(s1) 'vSub' num2str(s2) '   ',num2str(avg_withincluster_indistance_assigned_tasksimilaritymat(s1,s2))];  %save the data as a string to be written to the output
                dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
                
                
            end
        end
    end
end
            
            