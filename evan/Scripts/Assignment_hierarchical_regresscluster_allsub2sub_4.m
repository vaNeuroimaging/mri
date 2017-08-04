%subnetworks = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Consensusfillin_Powercolors_cleaned_rawassn_minsize400_regularize.dtseries.nii';
subnetworksfile = ['Template_match_systems_Poldrome_MSC.dtseries.nii'];
subparcelfile = ['Parcels_Poldrome_MSC.dtseries.nii'];



%'/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_LR_timeseries.dtseries.nii';
%dynamic_distance_metric_thresh = 1;
doublematch_factor = 10;
absolute_distance_cutoff = 1;%.85;
absolute_distance_cutoff_unregressed = .25;
xdistance = 20;

%parcelpct_inclusterthresh = .5;

MSCnames = {'MSC01','MSC02','MSC03','MSC04'};



% groupsurfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
% [ign, groupciftifiles] = textread(groupsurfdatafile,'%s %s');
% 
% grouptmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
% [groupsubjects, grouptmasks] = textread(grouptmaskfile,'%s %s');
% 
% groupsubstouse = [1:120];
% groupsubjects = groupsubjects(groupsubstouse);
% grouptmasks = grouptmasks(groupsubstouse);
% groupciftifiles = groupciftifiles(groupsubstouse);

%groupnetworks = ['/data/cn4/evan/RestingState/Ind_variability/Poldrome/Assignment/regresscluster_v' MSCname '/' MSCname '_Template_match.dtseries.nii'];%'/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii';
%groupparcelfile = ['/data/cn4/evan/RestingState/Ind_variability/Poldrome/Assignment/regresscluster_v' MSCname '/' MSCname '_parcels_LR_0.4.dtseries.nii'];%'/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii';
%groupdconn = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.dconn.nii';

ncortverts = 59412;
 


%distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
%distances = distances(1:ncortverts,1:ncortverts);

subnetworks = cifti_read(subnetworksfile); subnetworks = subnetworks(1:ncortverts,:);
subparcels = cifti_read(subparcelfile); subparcels = subparcels(1:ncortverts,:);
%subparcelIDs = unique(subparcels); subparcelIDs(subparcelIDs==0) = [];


subnetworks(subnetworks > 9) = 0; subnetworks(subnetworks < 9) = 0;
%subnetworks = subnetworks(:,1:2);
%subparcels = subparcels(:,1:2);

%groupparcels = cifti_read(groupparcelfile); groupparcels = groupparcels(1:ncortverts);
%groupparcelIDs = unique(groupparcels); groupparcelIDs(groupparcelIDs==0) = [];

clear adjacencies

parcels_subs_IDs = zeros(0,2);
for s = 1:size(subparcels,2);
    IDs = unique(subparcels(:,s));
    IDs(IDs==0) = [];
    parcels_subs_IDs = [parcels_subs_IDs ; [repmat(s,length(IDs),1) IDs]];
    
    adjacencies{s} = eye(length(IDs));
    
    borderverts = find(subparcels(:,s)==0);
    for vert = borderverts(:)'
        
        [vertdistances sorti] = sort(distances(vert,:),'ascend');
            
        sorti(subparcels(sorti,s)==0) = [];
        closestparcel = subparcels(sorti(1),s);
        closestparcelind = find(IDs==closestparcel);
            
        sorti(subparcels(sorti,s)==closestparcel) = [];
        secondclosestparcel = subparcels(sorti(1),s);
        secondclosestparcelind = find(IDs==secondclosestparcel);
        
        adjacencies{s}(closestparcelind,secondclosestparcelind) = 1;
        adjacencies{s}(secondclosestparcelind,closestparcelind) = 1;
    end
    
    
    
end
parcels_subs_IDs(:,end+1) = 0;

parcel_centroids = zeros(size(parcels_subs_IDs,1),1);

for i = 1:size(parcels_subs_IDs,1)
    thissub = parcels_subs_IDs(i,1);
    thisID = parcels_subs_IDs(i,2);
    parcelinds = find(subparcels(:,thissub)==thisID);
    networkIDs_byvert = subnetworks(parcelinds,thissub);
    networkIDs_byvert(networkIDs_byvert==0) = [];
    thisnetworkID = mode(networkIDs_byvert);
%    if (nnz(subnetworks(parcelinds,thissub)==thisnetworkID) / nnz(parcelinds)) > parcelpct_inclusterthresh
        parcels_subs_IDs(i,3) = thisnetworkID;
        
%         subparcels(parcelinds(subnetworks(parcelinds,thissub)~=thisnetworkID),thissub) = 0;
%         parcelinds = find(subparcels(:,thissub)==thisID);
        
%    end
    
    parcel_totaldistances = sum(distances(parcelinds,parcelinds),2);
    [ign mini] = min(parcel_totaldistances);
    parcel_centroids(i) = parcelinds(mini);
    
end

match_matrix = zeros(size(parcels_subs_IDs,1));


for s1 = 1:size(subparcels,2)
    
    if s1==1
        %    disp('Loading subject distances')
        %     distances = smartload('/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/normalwall_distmat_333/distmat_homo_surf_geodesic_vol_euc.mat');
        %     distances = distances(1:ncortverts,1:ncortverts);
        %     distances(1:29696,29697:end) = 1000; distances(29697:end,1:29696) = 1000;
        sub1timeseries = cifti_read('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_LR_timeseries.dtseries.nii');
        tmask = load('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_total_tmask.txt');
        sub1timeseries = sub1timeseries(1:ncortverts,logical(tmask));
    elseif s1 < size(subparcels,2)
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
    
    for s2 = 1:size(subparcels,2)
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
            
            
            
            
            disp(['Matching network clusters: subject ' num2str(s1) ' vs subject ' num2str(s2)'])
            [s1matches, s2matches] = networkClusters_match(subnetworks(:,s1),distances,subnetworks(:,s2));
            
            
            
            disp(['Calculating subject ' num2str(s1) ' parcel connectivity patterns regressing the system cluster'])
            
            
            sub1clustertimecourses = zeros(size(s1matches,2),size(sub1timeseries,2));
            for i = 1:size(s1matches,2)
                sub1clustertimecourses(i,:) = mean(sub1timeseries(logical(s1matches(:,i)),:),1);
            end
            
            sub1parcelIDs = parcels_subs_IDs(parcels_subs_IDs(:,1)==s1,2);
            sub1parcel_clusterconnections = zeros(length(sub1parcelIDs),size(s1matches,2));
            sub1parcel_unregressedclusterconnections = zeros(length(sub1parcelIDs),size(s1matches,2));
            
            sub1parcels_incluster = zeros(length(sub1parcelIDs),1);
            for i = 1:length(sub1parcelIDs)
            
                parceltimecourse = mean(sub1timeseries(subparcels(:,s1)==sub1parcelIDs(i),:),1);
                
                numparcelverts = nnz(subparcels(:,s1)==sub1parcelIDs(i));
                
                numparcelverts_incluster = zeros(size(s1matches,2),1);
                for clusternum = 1:size(s1matches,2)
                    numparcelverts_incluster(clusternum) = nnz((subparcels(:,s1)==sub1parcelIDs(i)) .* s1matches(:,clusternum));
                end
                [maxverts_inclus maxi] = max(numparcelverts_incluster);
                
                [modeval nummode] = mode(subnetworks(subparcels(:,s1)==sub1parcelIDs(i),s1));
                
                if maxverts_inclus == nummode
                    %(maxverts_inclus/numparcelverts) >= parcelpct_inclusterthresh
                    
                    sub1parcels_incluster(i) = maxi;
                    
                end
                if sub1parcels_incluster(i) > 0
                    sub1parcel_clusterconnections(i,:) = FisherTransform(partialcorr(parceltimecourse',sub1clustertimecourses',sub1clustertimecourses(sub1parcels_incluster(i),:)'));
                    sub1parcel_unregressedclusterconnections(i,:) = FisherTransform(paircorr_mod(parceltimecourse',sub1clustertimecourses'));
                end
            end
            
            
            disp(['Calculating subject ' num2str(s2) ' parcel connectivity patterns regressing the system cluster'])
            
            sub2clustertimecourses = zeros(size(s2matches,2),size(sub2timeseries,2));
            for i = 1:size(s2matches,2)
                sub2clustertimecourses(i,:) = mean(sub2timeseries(logical(s2matches(:,i)),:),1);
            end
            
            %subparcelpatterns = zeros(ncortverts,length(subparcelIDs));
            sub2parcelIDs = parcels_subs_IDs(parcels_subs_IDs(:,1)==s2,2);
            sub2parcel_clusterconnections = zeros(length(sub2parcelIDs),size(s2matches,2));
            sub2parcel_unregressedclusterconnections = zeros(length(sub2parcelIDs),size(s2matches,2));
            
            sub2parcels_incluster = zeros(length(sub2parcelIDs),1);
            for i = 1:length(sub2parcelIDs)
                
                parceltimecourse = mean(sub2timeseries(subparcels(:,s2)==sub2parcelIDs(i),:),1);
                
                numparcelverts = nnz(subparcels(:,s2)==sub2parcelIDs(i));
                
                numparcelverts_incluster = zeros(size(s2matches,2),1);
                for clusternum = 1:size(s2matches,2)
                    numparcelverts_incluster(clusternum) = nnz((subparcels(:,s2)==sub2parcelIDs(i)) .* s2matches(:,clusternum));
                end
                [maxverts_inclus maxi] = max(numparcelverts_incluster);
                
                [modeval nummode] = mode(subnetworks(subparcels(:,s2)==sub2parcelIDs(i),s2));
                
                if maxverts_inclus == nummode
                    %(maxverts_inclus/numparcelverts) >= parcelpct_inclusterthresh
                    
                    sub2parcels_incluster(i) = maxi;
                    
                end
                if sub2parcels_incluster(i) > 0
                    sub2parcel_clusterconnections(i,:) = FisherTransform(partialcorr(parceltimecourse',sub2clustertimecourses',sub2clustertimecourses(sub2parcels_incluster(i),:)'));
                    sub2parcel_unregressedclusterconnections(i,:) = FisherTransform(paircorr_mod(parceltimecourse',sub2clustertimecourses'));
                end
            end
            
            
            
            
            disp(['Assigning parcels'])
            
            matchnum = 0;
            
            s1inds_inmatchmatrix = find(parcels_subs_IDs(:,1)==s1);
            s2inds_inmatchmatrix = find(parcels_subs_IDs(:,1)==s2);
            
            for clusternum = 1%:size(s2matches,2);
                sub1parcelnums = find(sub1parcels_incluster==clusternum);
                sub2parcelnums = find(sub2parcels_incluster==clusternum);
                
                distance_metrics = zeros(length(sub1parcelnums),length(sub2parcelnums));
                if min(size(distance_metrics)) > 1
                    for i = 1:length(sub1parcelnums)
                        for j = 1:length(sub2parcelnums)
                            non_naninds = logical((~isnan(sub1parcel_clusterconnections(sub1parcelnums(i),:))) .* (~isnan(sub2parcel_clusterconnections(sub2parcelnums(j),:))));
                            distance_metrics(i,j) = 1-paircorr_mod(sub1parcel_clusterconnections(sub1parcelnums(i),non_naninds)',sub2parcel_clusterconnections(sub2parcelnums(j),non_naninds)');
                        end
                    end
                    
                    distance_metrics(distance_metrics > absolute_distance_cutoff) = Inf;
                    
                elseif min(size(distance_metrics)) == 1
                    for i = 1:length(sub1parcelnums)
                        for j = 1:length(sub2parcelnums)
                            non_naninds = logical((~isnan(sub1parcel_unregressedclusterconnections(sub1parcelnums(i),:))) .* (~isnan(sub2parcel_unregressedclusterconnections(sub2parcelnums(j),:))));
                            distance_metrics(i,j) = 1-paircorr_mod(sub1parcel_unregressedclusterconnections(sub1parcelnums(i),non_naninds)',sub2parcel_unregressedclusterconnections(sub2parcelnums(j),non_naninds)');
                        end
                    end
                    
                    distance_metrics(distance_metrics > absolute_distance_cutoff_unregressed) = Inf;
                end
                    
                if min(size(distance_metrics)) > 0
                    
                    geo_distmat = distances(parcel_centroids(s1inds_inmatchmatrix(sub1parcelnums)),parcel_centroids(s2inds_inmatchmatrix(sub2parcelnums)));
                    distance_metrics(geo_distmat > xdistance) = Inf;
                                
                    this_assign_matrix = zeros(size(distance_metrics));
                    
                    [ign mini] = min(distance_metrics,[],2);
                    s1matchinds = sub2ind(size(distance_metrics),[1:size(distance_metrics,1)]',mini);
                    [ign mini] = min(distance_metrics,[],1);
                    s2matchinds = sub2ind(size(distance_metrics),mini',[1:size(distance_metrics,2)]');
                    this_assign_matrix(union(s1matchinds,s2matchinds)) = 1;
                    this_assign_matrix(intersect(s1matchinds,s2matchinds)) = doublematch_factor;
                    this_assign_matrix(isinf(distance_metrics)) = 0;
                    
                    % eliminate non-adjacent sub1 parcel matches
                    for i = 1:size(this_assign_matrix,2)
                        parcels_that_matched_indswithinmat = find(this_assign_matrix(:,i));
                        matches_to_thisparcel = sub1parcelnums(parcels_that_matched_indswithinmat);
                        these_adjacencies = adjacencies{s1}(matches_to_thisparcel,matches_to_thisparcel);
                        if any(these_adjacencies(:)==0)
                        clusters = zeros(length(matches_to_thisparcel),1);
                            for matchnum = 1:length(matches_to_thisparcel)
                                adjacentmatches = find(these_adjacencies(matchnum,:));
                                adjacentclustervals = clusters(adjacentmatches);
                                adjacentclustervals(adjacentclustervals==0) = [];
                                if isempty(adjacentclustervals)
                                    clusters(adjacentmatches) = matchnum;
                                elseif length(adjacentclustervals)==1
                                    clusters(adjacentmatches) = adjacentclustervals;
                                else
                                    for value = 2:length(adjacentclustervals)
                                        clusters(clusters==adjacentclustervals(value)) = adjacentmatches(1);
                                    end
                                end
                            end
                            
                            geo_distances_toparcels = geo_distmat(parcels_that_matched_indswithinmat,i);
                            [ign mini] = min(geo_distances_toparcels);
                            %distances_toparcels = distance_metrics(parcels_that_matched_indswithinmat,i);
                            %[ign mini] = min(distances_toparcels);
                            closestcluster = clusters(mini);
                            this_assign_matrix(parcels_that_matched_indswithinmat(clusters~=closestcluster),i) = 0;
                        end
                    end
                    
                    % eliminate non-adjacent sub2 parcel matches
                    for i = 1:size(this_assign_matrix,1)
                        parcels_that_matched_indswithinmat = find(this_assign_matrix(i,:));
                        matches_to_thisparcel = sub2parcelnums(parcels_that_matched_indswithinmat);
                        these_adjacencies = adjacencies{s2}(matches_to_thisparcel,matches_to_thisparcel);
                        if any(these_adjacencies(:)==0)
                        clusters = zeros(length(matches_to_thisparcel),1);
                            for matchnum = 1:length(matches_to_thisparcel)
                                adjacentmatches = find(these_adjacencies(matchnum,:));
                                adjacentclustervals = clusters(adjacentmatches);
                                adjacentclustervals(adjacentclustervals==0) = [];
                                if isempty(adjacentclustervals)
                                    clusters(adjacentmatches) = matchnum;
                                elseif length(adjacentclustervals)==1
                                    clusters(adjacentmatches) = adjacentclustervals;
                                else
                                    for value = 2:length(adjacentclustervals)
                                        clusters(clusters==adjacentclustervals(value)) = adjacentmatches(1);
                                    end
                                end
                            end
                            
                            geo_distances_toparcels = geo_distmat(i,parcels_that_matched_indswithinmat);
                            [ign mini] = min(geo_distances_toparcels);
                            %distances_toparcels = distance_metrics(i,parcels_that_matched_indswithinmat);
                            %[ign mini] = min(distances_toparcels);
                            closestcluster = clusters(mini);
                            this_assign_matrix(i,parcels_that_matched_indswithinmat(clusters~=closestcluster)) = 0;
                        end
                    end
                    
                    
                            
                            
                     
                    
                    [this_s1matches this_s2matches] = find(this_assign_matrix);
                    for i = 1:length(this_s1matches)
                        this_s1ind = s1inds_inmatchmatrix(sub1parcelnums(this_s1matches(i)));
                        this_s2ind = s2inds_inmatchmatrix(sub2parcelnums(this_s2matches(i)));
                        
                        match_matrix(this_s1ind,this_s2ind) = this_assign_matrix(this_s1matches(i), this_s2matches(i));
                        match_matrix(this_s2ind,this_s1ind) = this_assign_matrix(this_s1matches(i), this_s2matches(i));
                    end
                    
                end
                
            end
        end
    end
end
    
mat2pajek_mod_EG(match_matrix,find(triu(match_matrix,1)),'pajek.net')
rawclrs = infomap_wrapper_mod_EG('pajek.net',100);
dlmwrite(['rawassn.txt'],rawclrs,'\t')
delete('pajek.clu')
simple_assigns = modify_clrfile('simplify','rawassn.txt',2);

indices_in_networks = find(parcels_subs_IDs(:,3) > 0);

[ign sorti] = sort(simple_assigns(indices_in_networks),'ascend');

sorted_inds_in_networks = indices_in_networks(sorti);

figure;
parcel_correlmat_figmaker(match_matrix(sorted_inds_in_networks,sorted_inds_in_networks),parcels_subs_IDs(sorted_inds_in_networks,3),[-2 2],'Match Matrix')

assigned_output= zeros(size(subparcels));

for ind = 1:length(simple_assigns)
    s = parcels_subs_IDs(ind,1);
    ID = parcels_subs_IDs(ind,2);
    assignment = simple_assigns(ind);
    assigned_output(subparcels(:,s)==ID,s) = assignment;
end

assigned_output(end+1:66697,:) = 0;

cifti_write_wHDR(assigned_output,[],'Assignments')
    