assignmentsfile = 'Assignments_xd15_1.0thr_merged.dtseries.nii';
%assignments = cifti_read(assignmentsfile);
assignmentsstruct = ft_read_cifti_mod(assignmentsfile); assignments = assignmentsstruct.data; assignmentsstruct.data = [];

IDs = unique(assignments); IDs(IDs<1) = [];

for ID = IDs(:)'
    disp(ID)
    
    IDclusters = zeros(66697,0);
    clustersubjects = [];
    multiples = 0;
    for s = 1:size(assignments,2)
        thisclusters = metric_cluster_cifti(assignments(:,s),ID-.01,ID+.01,0);
        IDclusters(:,end+1:end+size(thisclusters,2)) = thisclusters;
        clustersubjects(end+1:end+size(thisclusters,2)) = s;
        if size(thisclusters,2) > 1;
            multiples = multiples + 1;
        end
    end
    
    if multiples > 0
        centroids = zeros(size(IDclusters,2),1);
        for clusternum = 1:size(IDclusters,2)
            clusterinds = find(IDclusters(:,clusternum));
            [mindist mini] = min(sum(distances(clusterinds,clusterinds),2));
            centroids(clusternum) = clusterinds(mini);
        end
        
        centroid_distances = distances(centroids,centroids);
        
        if multiples == 1
            subject_withmultiples = mode(clustersubjects);
            multiple_log = logical(clustersubjects==subject_withmultiples);
            multiple_inds = find(multiple_log);
            
            total_distances_formultiples = sum(centroid_distances(multiple_log,~multiple_log),2);
            inds_todelete = multiple_inds(total_distances_formultiples > min(total_distances_formultiples));
            
            assignments(logical(sum(IDclusters(:,inds_todelete),2)),subject_withmultiples) = 0;
        
        else
            
            %samesub_mat = ~logical(squareform(pdist(clustersubjects'))) - eye(length(clustersubjects));
            %centroid_distances(logical(samesub_mat)) = 100;
            
            clustering = linkage(squareform(centroid_distances), 'average');
            clusters = cluster(clustering, 'MaxClust', [1:80]);

            for numclust = 1:size(clusters,2)
                Qvals(numclust) = M_calc_modularity(clusters(:,numclust),centroid_distances);
            end
            [maxQval maxQind] = max(Qvals);
            
            %cut = clustering(end-maxQind-1,3);
            
            
            %colormap = [1 0 0;0 0 1;1 1 0; 0 .8 0; 1 .6 1;1 .5 0;1 .7 .4;0 .6 .6;.6 .2 1];
            
            %[H,T,perm] = dendrogram_evan(clustering, colormap,1,0, 'orientation','left', 'colorthreshold', cut+.00001 ,'ORIENTATION','left');
            
            ordered_clusterIDs = clusters(:,maxQind+1);
            
            for ID = 2:max(ordered_clusterIDs);
                newassignmentval = max(unique(assignments)) + 1;
                clusterinds = find(ordered_clusterIDs==ID);
                for thisind = clusterinds(:)'
                    thissub = clustersubjects(thisind);
                    assignments(logical(IDclusters(:,thisind)),thissub) = newassignmentval;
                end
            end
            
        end
    end
end

outfilename = [assignmentsfile(1:end-13) '_adjacent'];
assignmentsstruct.data = assignments;
%cifti_write_wHDR(assignments,[],outfilename);
ft_write_cifti_mod(outfilename,assignmentsstruct);
        
        
        
        
        