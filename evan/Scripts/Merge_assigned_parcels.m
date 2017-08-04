ncortverts = 59412;

assignmentsname = 'Assignments_xd15_1.0thr.dtseries.nii';
assignments = cifti_read(assignmentsname); assignments = assignments(1:ncortverts,:);
parcels = cifti_read('Parcels_Poldrome_MSC.dtseries.nii'); parcels = parcels(1:ncortverts,:);
neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');



if ~exist('distances')
    distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
    distances = distances(1:ncortverts,1:ncortverts);
end

merged1 = assignments;

for s = 1:size(assignments,2)
    disp(['Subject ' num2str(s) ', first pass'])
    
    parcelIDs = unique(parcels(:,s)); parcelIDs(parcelIDs==0) = [];
    
    borderverts = find(parcels(:,s)==0);
    
    for vert = borderverts(:)'
        
        vertneighs = neighbors(vert,2:end); vertneighs(isnan(vertneighs)) = [];
        neighvals = assignments(vertneighs,s);
        
        if all(neighvals==mode(neighvals(neighvals>0)))
            merged1(vert,s) = neighvals(1);
        else
        
        neighvals(neighvals==0) = [];
        
        if numel(unique(neighvals)) < 2
            
            [vertdistances sorti] = sort(distances(vert,:),'ascend');
            
            sorti(parcels(sorti,s)==0) = [];
            closestparcel = parcels(sorti(1),s);
            
            sorti(parcels(sorti,s)==closestparcel) = [];
            secondclosestparcel = parcels(sorti(1),s);
            
            if mode(assignments(parcels(:,s)==closestparcel,s))==mode(assignments(parcels(:,s)==secondclosestparcel,s))
                
                merged1(vert,s) = mode(assignments(parcels(:,s)==closestparcel,s));
            end
        end
        end
    end
    
    merged2 = merged1;
    
    disp('second pass')
    
    tempclusters = metric_cluster_cifti(merged1(:,s),.001,max(merged1(:,s)) + .1,0);
    clusters = zeros(size(merged1(:,s)));
    for i = 1:size(tempclusters,2)
        clusters = clusters + tempclusters(:,i).*i;
    end
    unassignedparcels = parcels(:,s) .* (clusters==0); unassignedparcels(unassignedparcels>0) = unassignedparcels(unassignedparcels>0) + max(clusters);
    clusters = clusters + unassignedparcels;
    clusterIDs = unique(clusters); clusterIDs(clusterIDs==0) = [];
    
    borderverts = find(clusters==0);
    
    for vert = borderverts(:)'
        
        vertneighs = neighbors(vert,2:end); vertneighs(isnan(vertneighs)) = [];
        neighvals = merged1(vertneighs,s);
        
        neighvals(neighvals==0) = [];
        
        if numel(unique(neighvals)) < 2
            
            [vertdistances sorti] = sort(distances(vert,:),'ascend');
            
            sorti(clusters(sorti)==0) = [];
            closestparcel = clusters(sorti(1));
            
            sorti(clusters(sorti)==closestparcel) = [];
            secondclosestparcel = clusters(sorti(1));
            
            if mode(merged1(clusters==closestparcel,s))==mode(merged1(clusters==secondclosestparcel,s))
                
                merged2(vert,s) = mode(merged1(clusters==closestparcel,s));
            end
        end
        
    end
    
end

IDs = unique(merged2); IDs(IDs<1) = [];
recoloredIDs = randperm(max(IDs));
merged_recolored = zeros(size(merged2));
for IDnum = 1:length(IDs)
    merged_recolored(merged2 == IDs(IDnum)) = recoloredIDs(IDnum);
end

outname = [assignmentsname(1:end-13) '_merged.dtseries.nii'];
merged_recolored(end+1:66697,:) = 0;

cifti_write_wHDR(merged_recolored,[],outname);
            