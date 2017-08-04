function Cluster_and_combine(assignmentsfile,hem,minsubcluster)
%Cluster_and_combine(assignmentsfile,hem)

%hem = 'R';
%assignmentsfile = 'Subject_clustering_vertexwise_allnetworks_R_minQ_0.1.func.gii';%['Subject_clustering_vertexwise_' hem '.func.gii'];
sizesfile = [assignmentsfile(1:end-9) '_clustersize.func.gii'];
%minsubclusters = [10];


minclustersizemm = 20;

surfacearea = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR_surfaceareas.func.gii']);
surfacearea = surfacearea.cdata;


subclustersizes = gifti(sizesfile); subclustersizes = subclustersizes.cdata;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%for minsubcluster = minsubclusters
    assignments = gifti(assignmentsfile); assignments = assignments.cdata;
    
    assignments(subclustersizes < minsubcluster) = 0;
    
    assignments = assignments(:,1:3);
    
    networks = unique(assignments); networks(networks==0) = [];
    
    for networknum = 1:length(networks)
        
        networkpresent = any(assignments==networks(networknum),2);
        data_thisval = find(networkpresent);
        clusteredmetric = zeros(size(networkpresent));
        
        for vertex = data_thisval'
            
            %find the neighbors of this vertex
            vertexneighbors = neighbors(vertex,:);
            
            %find which of those neighbors also pass the thresholds
            vertexneighbors_inthresh = intersect(data_thisval,vertexneighbors);
            
            %find if those neighbors have already been assigned different cluster values
            uniqueneighborvals = unique(clusteredmetric(vertexneighbors_inthresh));
            uniqueneighborvals(uniqueneighborvals==0) = [];
            
            %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
            if isempty(uniqueneighborvals)
                clusteredmetric(vertexneighbors_inthresh) = vertex;
                %if there is only one previous cluster identifier present, make all the neighbors that value
            elseif length(uniqueneighborvals)==1
                clusteredmetric(vertexneighbors_inthresh) = uniqueneighborvals;
                %if there are multiple cluster identifier values in the neighborhood, merge them into one
            else
                for valuenum = 2:length(uniqueneighborvals)
                    clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
                end
            end
            
        end
        
        %find out what the unique cluster identifier values are
        uniqueclustervals = unique(clusteredmetric);
        uniqueclustervals(uniqueclustervals==0) = [];
        
        for clusternum = 1:length(uniqueclustervals)
            
            if sum(surfacearea(clusteredmetric==uniqueclustervals(clusternum))) < minclustersizemm
                %length(find(clusteredmetric==uniqueclustervals(clusternum))) < minclustersizemm
                
                indicestozero = find(clusteredmetric==uniqueclustervals(clusternum));
                
                for indexnum = indicestozero'
                    assignments(indexnum,(assignments(indexnum,:)==networks(networknum))) = 0;
                end
                
            end
        end
        
    end
    
    
    
    sorted_assignments = sort(assignments,2);
    [combinations ign uniqueIDs] = unique(sorted_assignments,'rows');
    
    
    save(gifti(single(uniqueIDs)),[assignmentsfile(1:end-9) '_clusters.func.gii'])
    %%
    system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample ' assignmentsfile(1:end-9) '_clusters.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  ' assignmentsfile(1:end-9) '_clusters_164.func.gii -largest'])
    
    combinations_upsampled = gifti([assignmentsfile(1:end-9) '_clusters_164.func.gii']); combinations_upsampled = combinations_upsampled.cdata;
    
    sphere = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii']);
    [phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    thetavals = -pi/2 : 1/70 : pi/2;
    surf_withlines = zeros(size(phi));
    
    for i = 1:length(thetavals)-1
        theseindices = intersect(find(theta>thetavals(i)), find(theta>thetavals(i+1)));
        surf_withlines(theseindices) = i;
    end
    
    surf_withlines = surf_withlines+1;
    
    %combinationIDs = unique(combinations_upsampled); combinationIDs(combinationIDs==0) = [];
    
    finaloutput = zeros(size(combinations_upsampled));
    listofuniqueIDs = unique(uniqueIDs);
    for combinationnum = 1:length(listofuniqueIDs)
        
        combination_withlines = surf_withlines .* (combinations_upsampled==listofuniqueIDs(combinationnum));
        lines_within_parcel = unique(combination_withlines);
        lines_within_parcel(lines_within_parcel==0)=[];
        
        if ~isempty(lines_within_parcel)
            
            values_toassign = sorted_assignments(uniqueIDs==listofuniqueIDs(combinationnum),:);
            values_toassign = values_toassign(1,:);
            values_toassign(values_toassign==0) = [];
            
            
            if length(values_toassign) <= length(lines_within_parcel)
                
                for comnum = 1:length(values_toassign)
                    
                    lineval_indices = [comnum : length(values_toassign) : length(lines_within_parcel)];
                    linevals_for_this_community = lines_within_parcel(lineval_indices);
                    for lineval = linevals_for_this_community'
                        finaloutput(combination_withlines==lineval) = values_toassign(comnum);
                    end
                    
                end
            end
        end
    end
    
    save(gifti(single(finaloutput)),[assignmentsfile(1:end-9) 'minsize_' num2str(minsubcluster) '_combineclusters_164.func.gii']);%['Subject_clustering_vertexwise_' hem '_combineclusters_164.func.gii']);
    
%end

