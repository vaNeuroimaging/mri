subnetworks = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Consensusfillin_Powercolors_cleaned_rawassn_minsize400_regularize.dtseries.nii';
subparcelfile = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_84_subsurf_edge_LR_watershedmerge_0.45_tweaked.dtseries.nii';
subtimeseriesfile = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_LR_timeseries.dtseries.nii';
subtmask = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_correlation_normalwall/84sub_333_all/allsubs_total_tmask.txt';
subdistmat = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/normalwall_distmat_333/distmat_homo_surf_geodesic_vol_euc.mat';

groupnetworks = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii';
groupparcelfile = '/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii';
groupdconn = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_normalwall/120_avg_corr_LR.dconn.nii';

groupsubstouse = [1:120];
ncortverts = 59412;
 
kden_thresh = .05;

xdistance = 30;

distance_metric_thresh = .5;

loadprevious = 1;

subparcels = cifti_read(subparcelfile); subparcels = subparcels(1:ncortverts);
subparcelIDs = unique(subparcels); subparcelIDs(subparcelIDs==0) = [];

groupparcels = cifti_read(groupparcelfile); groupparcels = groupparcels(1:ncortverts);
groupparcelIDs = unique(groupparcels); groupparcelIDs(groupparcelIDs==0) = [];

% surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
% [subjects, ciftifiles] = textread(surfdatafile,'%s %s');
% 
% tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
% [subjects, tmasks] = textread(tmaskfile,'%s %s');
% 
% subjects = subjects(groupsubstouse);%(81:100);%(121:end);
% tmasks = tmasks(groupsubstouse);%(81:100);%(121:end);
% ciftifiles = ciftifiles(groupsubstouse);%(81:100);%(121:end);

disp('Matching network clusters')
[groupmatches, submatches] = networkClusters_matchtogroup(subnetworks,[],groupnetworks);

if ~loadprevious

disp('Loading subject distances')
distances = smartload(subdistmat);
distances = distances(1:ncortverts,1:ncortverts);
distances(1:29696,29697:end) = 1000; distances(29697:end,1:29696) = 1000;

disp('Calculating subject parcel connectivity patterns')
subtimeseries = cifti_read(subtimeseriesfile);
tmask = load(subtmask);
subtimeseries = subtimeseries(1:ncortverts,logical(tmask));


%subparcelpatterns = zeros(ncortverts,length(subparcelIDs));
subparcel_clusterconnections = zeros(length(subparcelIDs),size(submatches,2));
for i = 1:length(subparcelIDs)
    disp(['Parcel ' num2str(i)])
    subparcelpattern = FisherTransform(paircorr_mod(subtimeseries',mean(subtimeseries(subparcels==subparcelIDs(i),:),1)'));
    subparcelpattern(isnan(subparcelpattern)) = 0;
    distantinds = all((distances(:,subparcels==subparcelIDs(i)) > xdistance),2);
    for j = 1:size(submatches,2)
        indstoaverage = logical(distantinds .* submatches(:,j));
        subparcel_clusterconnections(i,j) = mean(subparcelpattern(indstoaverage));
    end
end
clear subtimeseries distances

disp('Loading group distances')
distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
distances = distances(1:ncortverts,1:ncortverts);
distances(1:29696,29697:end) = 1000; distances(29697:end,1:29696) = 1000;

disp('Loading group dense connectivity')

groupcorrel = cifti_read(groupdconn);
groupcorrel = groupcorrel(1:ncortverts,1:ncortverts);
groupcorrel(isnan(groupcorrel)) = 0;
%groupparcelpatterns = zeros(ncortverts,length(groupparcelIDs));

disp('Calculating group parcel connectivity patterns')
groupparcel_clusterconnections = zeros(length(groupparcelIDs),size(groupmatches,2));
for i = 1:length(groupparcelIDs)
    disp(['Parcel ' num2str(i)])
    groupparcelpattern = mean(groupcorrel(:,groupparcels==groupparcelIDs(i)),2);
    distantinds = all((distances(:,groupparcels==groupparcelIDs(i)) > xdistance),2);
    for j = 1:size(groupmatches,2)
        indstoaverage = logical(distantinds .* groupmatches(:,j));
        groupparcel_clusterconnections(i,j) = mean(groupparcelpattern(indstoaverage));
    end
end

clear groupcorrel distances

save('Group_sub_parcelcorrelpatterns.mat','groupparcel_clusterconnections','subparcel_clusterconnections','-v7.3')

else
    load('Group_sub_parcelcorrelpatterns.mat')
end
    

disp('Assigning parcels')

matchnum = 0;
submatchedparcels = zeros(ncortverts,1);
groupmatchedparcels = zeros(ncortverts,1);

for clusternum = 39%1:size(groupmatches,2);
    disp(['Network cluster number ' num2str(clusternum)])
    groupparcelnums = [];
    for i = 1:length(groupparcelIDs)
        numparcelverts = nnz(groupparcels==groupparcelIDs(i));
        numparcelverts_incluster = nnz((groupparcels==groupparcelIDs(i)) .* groupmatches(:,clusternum));
        if (numparcelverts_incluster / numparcelverts) > .5
            groupparcelnums = [groupparcelnums i];
        end
    end
    
    subparcelnums = [];
    for i = 1:length(subparcelIDs)
        numparcelverts = nnz(subparcels==subparcelIDs(i));
        numparcelverts_incluster = nnz((subparcels==subparcelIDs(i)) .* submatches(:,clusternum));
        if (numparcelverts_incluster / numparcelverts) > .5
            subparcelnums = [subparcelnums i];
        end
    end
    
    %distance_metrics = pdist(groupparcel_clusterconnections(groupparcelnums,:)',subparcel_clusterconnections(subparcelnums,:)');
    %distance_metrics = 1-(paircorr_mod(groupparcel_clusterconnections(groupparcelnums,:)',subparcel_clusterconnections(subparcelnums,:)'));
    distance_metrics = zeros(length(groupparcelnums),length(subparcelnums));
    for i = 1:length(groupparcelnums)
        for j = 1:length(subparcelnums)
            non_naninds = logical((~isnan(groupparcel_clusterconnections(groupparcelnums(i),:))) .* (~isnan(subparcel_clusterconnections(subparcelnums(j),:))));
            
            distance_metrics(i,j) = 1-paircorr_mod(groupparcel_clusterconnections(groupparcelnums(i),non_naninds)',subparcel_clusterconnections(subparcelnums(j),non_naninds)');
        end
    end
    
    distance_metrics(distance_metrics > distance_metric_thresh) = Inf;
    
    [subtogroup_assign,grouptosub_assign,subtogroup_cost, grouptosub_cost] = munkres_mult(distance_metrics,5);
    
    unique_subtogroup = unique(subtogroup_assign); unique_subtogroup(unique_subtogroup==0) = [];
    thiscluster_groupmatchparcels = zeros(ncortverts,length(unique_subtogroup));
    thiscluster_submatchparcels = zeros(ncortverts,length(unique_subtogroup));
    for i = 1:length(unique_subtogroup)
        groupparcels_inthiscluster{i} = find(subtogroup_assign==unique_subtogroup(i));
        for j = groupparcels_inthiscluster{i}(:)';
            thiscluster_groupmatchparcels(:,i) = thiscluster_groupmatchparcels(:,i) + (groupparcels==groupparcelIDs(groupparcelnums(j)));
        end
    end
    unique_grouptosub = unique(grouptosub_assign); unique_grouptosub(unique_grouptosub==0) = [];
    for j = 1:length(unique_grouptosub)
        for i = 1:length(groupparcels_inthiscluster)
            if any(groupparcels_inthiscluster{i}==unique_grouptosub(j))
                targetcluster = i;
                break
            end
        end
        subparcels_inthiscluster{j} = find(grouptosub_assign==unique_grouptosub(j));
        for k = subparcels_inthiscluster{j}(:)'
            thiscluster_submatchparcels(:,targetcluster) = thiscluster_submatchparcels(:,targetcluster) + (subparcels==subparcelIDs(subparcelnums(k)));
        end
           
    end
    clear groupclusters_inthiscluster subparcels_inthiscluster
    
    inds_toremove = logical(all(thiscluster_submatchparcels==0,1) + all(thiscluster_groupmatchparcels==0,1));
    thiscluster_groupmatchparcels(:,inds_toremove) = [];
    thiscluster_submatchparcels(:,inds_toremove) = [];
    
    for i = 1:size(thiscluster_submatchparcels,2)
        matchnum = matchnum+1;
        groupmatchedparcels(logical(thiscluster_groupmatchparcels(:,i))) = matchnum;
        submatchedparcels(logical(thiscluster_submatchparcels(:,i))) = matchnum;
    end
end

groupmatchedparcels(end+1:66697,:) = 0;
cifti_write_wHDR(groupmatchedparcels,[],'Groupmatchedparcels')
submatchedparcels(end+1:66697,:) = 0;
cifti_write_wHDR(submatchedparcels,[],'Submatchedparcels')
    
    %groupmatches(:,(end+1):(end+size(thiscluster_groupmatchparcels,2))) = thiscluster_groupmatchparcels;
    %submatches(:,(end+1):(end+size(thiscluster_submatchparcels,2))) = thiscluster_submatchparcels;
    