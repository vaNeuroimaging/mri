%subnetworks = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Consensusfillin_Powercolors_cleaned_rawassn_minsize400_regularize.dtseries.nii';
subnetworksfile = ['Template_match_systems_Poldrome_MSC.dtseries.nii'];
subparcelfile = ['Parcels_Poldrome_MSC.dtseries.nii'];



%'/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_LR_timeseries.dtseries.nii';
%dynamic_distance_metric_thresh = 1;
doublematch_factor = 10;
absolute_distance_cutoff = .5;
dicethresh = .25;

parcelpct_inclusterthresh = .5;

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
 
xdistance = 30;

%distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
%distances = distances(1:ncortverts,1:ncortverts);

subnetworks = cifti_read(subnetworksfile); subnetworks = subnetworks(1:ncortverts,:);
subparcels = cifti_read(subparcelfile); subparcels = subparcels(1:ncortverts,:);
%subparcelIDs = unique(subparcels); subparcelIDs(subparcelIDs==0) = [];


subnetworks(subnetworks > 3) = 0; subnetworks(subnetworks < 3) = 0;
%subnetworks = subnetworks(:,1:2);
%subparcels = subparcels(:,1:2);

%groupparcels = cifti_read(groupparcelfile); groupparcels = groupparcels(1:ncortverts);
%groupparcelIDs = unique(groupparcels); groupparcelIDs(groupparcelIDs==0) = [];


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
    networkIDs_byvert = subnetworks(parcelinds,thissub);
    networkIDs_byvert(networkIDs_byvert==0) = [];
    thisnetworkID = mode(networkIDs_byvert);
    if (nnz(subnetworks(parcelinds,thissub)==thisnetworkID) / nnz(parcelinds)) > parcelpct_inclusterthresh
        parcels_subs_IDs(i,3) = thisnetworkID;
    end
    
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
    else
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
            
            
            %subclusterpatterns = zeros(ncortverts,size(submatches,2));
            sub1clustertimecourses = zeros(size(s1matches,2),size(sub1timeseries,2));
            for i = 1:size(s1matches,2)
                %subclusterpatterns(:,i) = FisherTransform(paircorr_mod(subtimeseries',mean(subtimeseries(logical(submatches(:,i)),:),1)'));
                sub1clustertimecourses(i,:) = mean(sub1timeseries(logical(s1matches(:,i)),:),1);
            end
            %subclusterpatterns(isnan(subclusterpatterns)) = 0;
            
            %subparcelpatterns = zeros(ncortverts,length(subparcelIDs));
            sub1parcelIDs = parcels_subs_IDs(parcels_subs_IDs(:,1)==s1,2);
            sub1parcel_clusterconnections = zeros(length(sub1parcelIDs),size(s1matches,2));
            
            sub1parcels_incluster = zeros(length(sub1parcelIDs),1);
            regressed_parceltimecourses = rand([size(sub1timeseries,2),length(sub1parcelIDs)]);
            for i = 1:length(sub1parcelIDs)
                %disp(['Parcel ' num2str(i)])
                
                parceltimecourse = mean(sub1timeseries(subparcels(:,s1)==sub1parcelIDs(i),:),1);
                
                numparcelverts = nnz(subparcels(:,s1)==sub1parcelIDs(i));
                
                numparcelverts_incluster = zeros(size(s1matches,2),1);
                for clusternum = 1:size(s1matches,2)
                    numparcelverts_incluster(clusternum) = nnz((subparcels(:,s1)==sub1parcelIDs(i)) .* s1matches(:,clusternum));
                end
                [maxverts maxi] = max(numparcelverts_incluster);
                
                if (maxverts/numparcelverts) >= parcelpct_inclusterthresh
                    
                    sub1parcels_incluster(i) = maxi;
                    
                    [ign,ign2,regressed_parceltimecourses(:,i)] = regress(parceltimecourse',sub1clustertimecourses(sub1parcels_incluster(i),:)');
                end
            end
            
            sub1regressedparcelpatterns = FisherTransform(paircorr_mod(sub1timeseries',regressed_parceltimecourses));
            sub1regressedparcelpatterns(isnan(sub1regressedparcelpatterns)) = 0;
            for i = 1:length(sub1parcelIDs)
                distantinds = all((distances(:,subparcels(:,s1)==sub1parcelIDs(i)) > xdistance),2);
                for j = 1:size(s1matches,2)
                    indstoaverage = logical(distantinds .* s1matches(:,j));
                    sub1parcel_clusterconnections(i,j) = mean(sub1regressedparcelpatterns(indstoaverage,i));
                end
                sub1regressedparcelpatterns(subparcels(:,s1)==sub1parcelIDs(i),i) = 0;
            end
            sorted = sort(sub1parcel_clusterconnections(:),'descend');
            sub1regressedparcelpatterns_binarized = (sub1regressedparcelpatterns >= sorted(round(length(sorted)*dicethresh)));
            
            
            
            %         subregressedparcelpatterns = FisherTransform(paircorr_mod(subtimeseries',regressed_parceltimecourse));
            %
            %         distantinds = all((distances(:,subparcels==subparcelIDs(i)) > xdistance),2);
            %         for j = 1:size(submatches,2)
            %             indstoaverage = logical(distantinds .* submatches(:,j));
            %             subparcel_clusterconnections(i,j) = mean(subregressedparcelpatterns(indstoaverage));
            %         end
            %         subregressedparcelpatterns(subparcels==subparcelIDs(i)) = 0;
            %         subregressed_allparcelpatterns(1:ncortverts,i) = subregressedparcelpatterns(1:ncortverts);
            
            % out = zeros(66697,size(sub1regressedparcelpatterns,2));
            % out(1:ncortverts,:) = sub1regressedparcelpatterns;
            %
            % cifti_write_wHDR(out,[],'Subregressed_allparcelpatterns')
            
            %clear sub1regressed_allparcelpatterns
            
            % disp('Loading group distances')
            % distances = smartload(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/normalwall_distmat_333/distmat_surf_geodesic_vol_euc.mat']);
            % distances = distances(1:ncortverts,1:ncortverts);
            % distances(1:29696,29697:end) = 1000; distances(29697:end,1:29696) = 1000;
            
            disp(['Calculating subject ' num2str(s2) ' parcel connectivity patterns regressing the system cluster'])
            
            %subclusterpatterns = zeros(ncortverts,size(submatches,2));
            sub2clustertimecourses = zeros(size(s2matches,2),size(sub2timeseries,2));
            for i = 1:size(s2matches,2)
                %subclusterpatterns(:,i) = FisherTransform(paircorr_mod(subtimeseries',mean(subtimeseries(logical(submatches(:,i)),:),1)'));
                sub2clustertimecourses(i,:) = mean(sub2timeseries(logical(s2matches(:,i)),:),1);
            end
            %subclusterpatterns(isnan(subclusterpatterns)) = 0;
            
            %subparcelpatterns = zeros(ncortverts,length(subparcelIDs));
            sub2parcelIDs = parcels_subs_IDs(parcels_subs_IDs(:,1)==s2,2);
            sub2parcel_clusterconnections = zeros(length(sub2parcelIDs),size(s2matches,2));
            
            sub2parcels_incluster = zeros(length(sub2parcelIDs),1);
            regressed_parceltimecourses = rand([size(sub2timeseries,2),length(sub2parcelIDs)]);
            for i = 1:length(sub2parcelIDs)
                %disp(['Parcel ' num2str(i)])
                
                parceltimecourse = mean(sub2timeseries(subparcels(:,s2)==sub2parcelIDs(i),:),1);
                
                numparcelverts = nnz(subparcels(:,s2)==sub2parcelIDs(i));
                
                numparcelverts_incluster = zeros(size(s2matches,2),1);
                for clusternum = 1:size(s2matches,2)
                    numparcelverts_incluster(clusternum) = nnz((subparcels(:,s2)==sub2parcelIDs(i)) .* s2matches(:,clusternum));
                end
                [maxverts maxi] = max(numparcelverts_incluster);
                
                if (maxverts/numparcelverts) >= parcelpct_inclusterthresh
                    
                    sub2parcels_incluster(i) = maxi;
                    
                    [ign,ign2,regressed_parceltimecourses(:,i)] = regress(parceltimecourse',sub2clustertimecourses(sub2parcels_incluster(i),:)');
                end
            end
            
            sub2regressedparcelpatterns = FisherTransform(paircorr_mod(sub2timeseries',regressed_parceltimecourses));
            sub2regressedparcelpatterns(isnan(sub2regressedparcelpatterns)) = 0;
            for i = 1:length(sub2parcelIDs)
                distantinds = all((distances(:,subparcels(:,s2)==sub2parcelIDs(i)) > xdistance),2);
                for j = 1:size(s2matches,2)
                    indstoaverage = logical(distantinds .* s2matches(:,j));
                    sub2parcel_clusterconnections(i,j) = mean(sub2regressedparcelpatterns(indstoaverage,i));
                end
                sub2regressedparcelpatterns(subparcels(:,s2)==sub2parcelIDs(i),i) = 0;
            end
            sorted = sort(sub2parcel_clusterconnections(:),'descend');
            sub2regressedparcelpatterns_binarized = (sub2regressedparcelpatterns >= sorted(round(length(sorted)*dicethresh)));
            
            % out = zeros(66697,size(groupregressedparcelpatterns,2));
            % out(1:ncortverts,:) = groupregressedparcelpatterns;
            %
            % cifti_write_wHDR(out,[],'Groupregressed_allparcelpatterns')
            
            % clear groupcorrel distances groupregressed_allparcelpatterns out
            %
            % save('Group_sub_parcelcorrelpatterns.mat','groupparcel_clusterconnections','subparcel_clusterconnections','groupparcels_incluster','subparcels_incluster','-v7.3')
            %
            % else
            %     load('Group_sub_parcelcorrelpatterns.mat')
            % end
            
            
            
            disp(['Assigning parcels'])
            
            matchnum = 0;
            %sub1matchedparcels = zeros(ncortverts,1);
            %sub2matchedparcels = zeros(ncortverts,1);
            
            s1inds_inmatchmatrix = find(parcels_subs_IDs(:,1)==s1);
            s2inds_inmatchmatrix = find(parcels_subs_IDs(:,1)==s2);
            
            for clusternum = 1:size(s2matches,2);
                %disp(['Network cluster number ' num2str(clusternum)])
                sub1parcelnums = find(sub1parcels_incluster==clusternum);
                sub2parcelnums = find(sub2parcels_incluster==clusternum);
                
                %distance_metrics = pdist(groupparcel_clusterconnections(groupparcelnums,:)',subparcel_clusterconnections(subparcelnums,:)');
                %distance_metrics = 1-(paircorr_mod(groupparcel_clusterconnections(groupparcelnums,:)',subparcel_clusterconnections(subparcelnums,:)'));
                distance_metrics = zeros(length(sub1parcelnums),length(sub2parcelnums));
                if min(size(distance_metrics)) > 1
                    for i = 1:length(sub1parcelnums)
                        for j = 1:length(sub2parcelnums)
                            non_naninds = logical((~isnan(sub1parcel_clusterconnections(sub1parcelnums(i),:))) .* (~isnan(sub2parcel_clusterconnections(sub2parcelnums(j),:))));
                            
                            s1parcel_clusterconnections = sub1regressedparcelpatterns_binarized(sub1parcelnums(i),non_naninds);
                            s2parcel_clusterconnections = sub2regressedparcelpatterns_binarized(sub2parcelnums(j),non_naninds);
                            distance_metrics(i,j) = 1-(nnz(s1parcel_clusterconnections .* s2parcel_clusterconnections) ./ nnz(s1parcel_clusterconnections + s2parcel_clusterconnections));
                            %1-paircorr_mod(sub1parcel_clusterconnections(sub1parcelnums(i),non_naninds)',sub2parcel_clusterconnections(sub2parcelnums(j),non_naninds)');
                        end
                    end
                    %dynamic_distance_metric_thresh = max([min(distance_metrics,[],1) min(distance_metrics,[],2)']);
                    %distance_metrics(distance_metrics > dynamic_distance_metric_thresh) = Inf;
                
                    
                    geo_distmat = distances(parcel_centroids(s1inds_inmatchmatrix(sub1parcelnums)),parcel_centroids(s2inds_inmatchmatrix(sub2parcelnums)));
                    distance_metrics(geo_distmat > xdistance) = Inf;
                    distance_metrics(distance_metrics > absolute_distance_cutoff) = Inf;
                    
                    this_assign_matrix = zeros(size(distance_metrics));
                    
                    [ign mini] = min(distance_metrics,[],2);
                    s1matchinds = sub2ind(size(distance_metrics),[1:size(distance_metrics,1)]',mini);
                    [ign mini] = min(distance_metrics,[],1);
                    s2matchinds = sub2ind(size(distance_metrics),mini',[1:size(distance_metrics,2)]');
                    this_assign_matrix(union(s1matchinds,s2matchinds)) = 1;
                    this_assign_matrix(intersect(s1matchinds,s2matchinds)) = doublematch_factor;
                    this_assign_matrix(isinf(distance_metrics)) = 0;
                    
                else
                    
                    this_assign_matrix = ones(size(distance_metrics));
                end
                    
                    %thresh = max(distance_metrics(logical(this_assign_matrix))) * dynamic_distance_metric_thresh;
                    %this_assign_matrix(distance_metrics <= thresh) = 1;
                    
%                     for i = 1:size(this_assign_matrix,1)
%                         if any(this_assign_matrix(i,:))
%                             thresh = min(distance_metrics(i,logical(this_assign_matrix(i,:)))) ./ dynamic_distance_metric_thresh;
%                             this_assign_matrix(i,distance_metrics(i,:) < thresh) = 1;
%                         end
%                     end
%                     
%                     for i = 1:size(this_assign_matrix,2)
%                         if any(this_assign_matrix(:,i))
%                             thresh = min(distance_metrics(logical(this_assign_matrix(:,i)),i)) ./ dynamic_distance_metric_thresh;
%                             this_assign_matrix(distance_metrics(:,i) < thresh,i) = 1;
%                         end
%                     end
                    
                    [this_s1matches this_s2matches] = find(this_assign_matrix);
                    for i = 1:length(this_s1matches)
                        this_s1ind = s1inds_inmatchmatrix(sub1parcelnums(this_s1matches(i)));
                        this_s2ind = s2inds_inmatchmatrix(sub2parcelnums(this_s2matches(i)));
                        
                        match_matrix(this_s1ind,this_s2ind) = this_assign_matrix(this_s1matches(i), this_s2matches(i));
                        match_matrix(this_s2ind,this_s1ind) = this_assign_matrix(this_s1matches(i), this_s2matches(i));
                    end
                    
                
                
                
                %[s2assigns , s1assigns] = ManytoManyAssignment(distance_metrics);
%                 [s2assigns , s1assigns,subtogroup_cost, grouptosub_cost] = munkres_mult(distance_metrics,7);
%                 
%                 for i = 1:length(s1assigns)
%                     if s1assigns(i) > 0
%                         this_s1ind = s1inds_inmatchmatrix(sub1parcelnums(i));
%                         this_s2ind = s2inds_inmatchmatrix(sub2parcelnums(s1assigns(i)));
%                         match_matrix(this_s1ind,this_s2ind) = 1;
%                         match_matrix(this_s2ind,this_s1ind) = 1;
%                     end
%                 end
%                 for i = 1:length(s2assigns)
%                     if s2assigns(i) > 0
%                         this_s1ind = s1inds_inmatchmatrix(sub1parcelnums(s2assigns(i)));
%                         this_s2ind = s2inds_inmatchmatrix(sub2parcelnums(i));
%                         match_matrix(this_s1ind,this_s2ind) = 1;
%                         match_matrix(this_s2ind,this_s1ind) = 1;
%                     end
%                 end
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
    