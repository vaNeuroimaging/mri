function Consensus_networkIDs_fillin(regularizedfilename,distances)

%Cost that determines a bad match
costcutoff = 20;

overlapcutoff = 30;

%Smallest acceptable contiguous cluster
sizethresh = 15;

%Template
groupnetworksfile = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR.dtseries.nii';

%Find file name
slashlocs = strfind(regularizedfilename,'/');
if isempty(slashlocs)
    regularizednetworksfolder = pwd;
    regularizednetworksfile = regularizedfilename;
else
    regularizednetworksfolder = regularizedfilename(1:slashlocs(end));
    regularizednetworksfile = regularizedfilename(slashlocs(end)+1 : end);
end


nsurfverts = 29696 + 29716;

%Load distances file if needed
if ~exist('distances')
    load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
    distances = distances(1:nsurfverts,1:nsurfverts);
end

%Colors to use for extra unmatched networks
extraassignments = [11.5 4 6 17 2.5 4.5 18:100];

%Load subject networks
if strcmp(regularizednetworksfile(end-3:end),'.txt')
    subnetworks = load(regularizednetworksfile);
    regularizednetworksfile = regularizednetworksfile(1:end-4);
else
    subnetworks = cifti_read(regularizednetworksfile);
end
subnetworks = subnetworks(1:nsurfverts,:);

%Deal with regularization by Jonathan's code
%if min(subnetworks(:)) == 1;
    subnetworks = subnetworks - 1;
%end

%Read group template
groupnetworks = cifti_read(groupnetworksfile);
groupnetworks(groupnetworks==13) = 0; groupnetworks(groupnetworks==14) = 0;
surfgroupnetworks = groupnetworks(1:nsurfverts);
groupIDs = unique(surfgroupnetworks); groupIDs(groupIDs<1) = [];
for i = 1:length(groupIDs)
    vertsinthisgroupID(i) = nnz(surfgroupnetworks);
end


filledsubnetworks = zeros(size(subnetworks));


maxnetworksmatched = 0;

%Loop through columns
for col = 1:size(subnetworks,2)
    disp(['Column ' num2str(col) ' of ' num2str(size(subnetworks,2))])
    surfsubnetworks = subnetworks(:,col);
    
    %Fill in unassigned verts with info from the next column up
    unassigned = find(surfsubnetworks<1);
    for unassignedindex = unassigned'
        thisassignments = subnetworks(unassignedindex,col:end);
        thisassignments(thisassignments<1) = [];
        if ~isempty(thisassignments)
            surfsubnetworks(unassignedindex) = thisassignments(1);
        end
    end
    
    filledsubnetworks(:,col) = surfsubnetworks;
    
    %Get IDs of subject networks at this column
    subIDs = unique(surfsubnetworks); subIDs(subIDs<1) = [];
    
    %Compare each group network to each subject network and calculate the
    %cost for a match
    mean_mindists = zeros(length(groupIDs),length(subIDs));
    noverlapverts = zeros(length(groupIDs),length(subIDs));
    for groupIDnum = 1:length(groupIDs)
        groupID = groupIDs(groupIDnum);
        
        for subIDnum = 1:length(subIDs)
            subID = subIDs(subIDnum);
            submindists = min(distances(surfgroupnetworks==groupID,surfsubnetworks==subID),[],1);
            submindists(submindists > 100) = [];
            groupmindists = min(distances(surfgroupnetworks==groupID,surfsubnetworks==subID),[],2);
            mean_mindists(groupIDnum,subIDnum) = mean([submindists(:); groupmindists(:)]);
            
            noverlapverts(groupIDnum,subIDnum) = nnz((surfgroupnetworks==groupID) .* (surfsubnetworks==subID));
            
        end
        
    end
    
    %Prevent high-cost matches 
    %mean_mindists(mean_mindists>costcutoff) = Inf;
    
    %mean_mindists(noverlapverts>overlapcutoff) = Inf;
    
    %Assign via Hungarian algorithm
    [colassigns{col},ign] = munkres(mean_mindists);
    for i = 1:length(colassigns{col})
        if noverlapverts(i,colassigns{col}(i)) < overlapcutoff
            colassigns{col}(i) = 0;
        end
    end
    
    allsubIDs{col} = subIDs;
    unassignedIDs = allsubIDs{col};
    all_mean_mindists{col} = mean_mindists;
    
    %if a group network doesn't get assigned, make the cost = the cost
    %cutoff
%     assignedcosts = zeros(size(colassigns{col}));
%     for i = 1:length(assignedcosts)
%         
%         if (colassigns{col}(i)) > 0 && (mean_mindists(i,colassigns{col}(i)) < costcutoff)
%             assignedcosts(i) = mean_mindists(i,colassigns{col}(i));
%             unassignedIDs(unassignedIDs==allsubIDs{col}(colassigns{col}(i))) = [];
%         else
%             assignedcosts(i) = costcutoff;
%         end
%         
%     end
    
    
    %Display some things
%     disp(['Assignments: ' num2str(colassigns{col})])
%     
%     disp(' ')
%     
%     if nnz(colassigns{col}) > maxnetworksmatched
%         maxnetworksmatched = nnz(colassigns{col});
%         bestmatch_column = col;
%     elseif (nnz(colassigns{col}) == maxnetworksmatched) && (sum(vertsinthisgroupID(logical(colassigns{col}))) >= sum(vertsinthisgroupID(logical(colassigns{bestmatch_column}))))
%         bestmatch_column = col;
%     end
    
    
    
end

allmatched = zeros(size(subnetworks));
matched = zeros(size(subnetworks,1),1);

for col = 1:size(subnetworks,2)
    
    %Map the minimum cost color assignments onto the brain
    assign = colassigns{col};
    unassignedIDs = allsubIDs{col};
    
    for assignednum = 1:length(assign)
        if (assign(assignednum)) > 0 && (all_mean_mindists{col}(assignednum,assign(assignednum)) < costcutoff)
            matched(filledsubnetworks(:,col)==allsubIDs{col}(assign(assignednum))) = groupIDs(assignednum);
            unassignedIDs(unassignedIDs==allsubIDs{col}(assign(assignednum))) = [];
        end
    end
    
    
    
    vertstofillin = [];
    for extraIDnum = 1:length(unassignedIDs)
        vertstofillin = [vertstofillin; find(filledsubnetworks(:,col)==unassignedIDs(extraIDnum))];
    end
    
    subIDsassigned = zeros(length(assign),1);
    for i = 1:length(assign)
        if assign(i) > 0
            subIDsassigned(i) = allsubIDs{col}(assign(i));
        end
    end
    
    assign_nozeros = assign; assign_nozeros(assign_nozeros==0) = [];
    stillunassignedIDs = [];
    vertstoreplace = zeros(size(matched));
    for vert = vertstofillin'
        for thiscol = col+1:size(filledsubnetworks,2)
            if any(filledsubnetworks(vert,thiscol)==allsubIDs{col}(assign_nozeros))
                matched(vert) = groupIDs(subIDsassigned==filledsubnetworks(vert,thiscol));
                break
            end
            if thiscol==size(filledsubnetworks,2)
                stillunassignedIDs = unique([stillunassignedIDs filledsubnetworks(vert,col)]);
                vertstoreplace(vert) = 1;
            end
        end
    end
    
    %Map unmatched networks onto the brain
    recoloredverts = zeros(size(matched));
    for extraIDnum = 1:length(stillunassignedIDs)
        recoloredverts(filledsubnetworks(:,col)==stillunassignedIDs(extraIDnum)) = extraassignments(extraIDnum);
    end
    matched(logical(vertstoreplace)) = recoloredverts(logical(vertstoreplace));
    
    allmatched(:,col) = matched;
    
end


%Get minimum cost column
% matchmetrics = zeros(size(allmatched,2),1);
% for col = 1:size(allmatched,2)
%     vertdists = [];
%     for groupIDnum = 1:length(groupIDs)
%         groupID = groupIDs(groupIDnum);
%         
%         submindists = min(distances(surfgroupnetworks==groupID,allmatched(:,col)==groupID),[],1);
%         submindists(submindists > 100) = [];
%         groupmindists = min(distances(surfgroupnetworks==groupID,allmatched(:,col)==groupID),[],2);
%         
%         vertdists = [vertdists; groupmindists; submindists(:)];
%     end
%     matchmetrics(col) = mean(vertdists);
% end
% 
% [ign, bestmatch_column] = min(matchmetrics);

matchmetrics = zeros(size(allmatched,2),1);
for col = 1:size(allmatched,2)
    matchmetrics(col) = nnz(allmatched(:,col)==surfgroupnetworks);
end
[ign, bestmatch_column] = max(matchmetrics);

bestmatch = allmatched(:,bestmatch_column);

disp(['Best matching column is ' num2str(bestmatch_column)])



%Clean up best match and write it out
cleaned = zeros(size(bestmatch));
hems = {'L','R'};
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']; maskL = gifti(maskname); maskL = ~maskL.cdata;
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    
    hemdata = zeros(32492,1);
    hemdata(logical(mask)) = bestmatch([1:nnz(mask)] + (nnz(maskL) * (hemnum-1)));
    
    allcolors = unique(hemdata(hemdata > 0));
    
    for color = allcolors'
        clusteredmetric = zeros(size(hemdata));
        thiscolorverts = find(hemdata==color);
        
        for vertex = thiscolorverts'
            
            %find the neighbors of this vertex
            vertexneighbors = neighbors(vertex,:);
            
            %find which of those neighbors also pass the thresholds
            vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
            
            %find if those neighbors have already been assigned different cluster values
            uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
            uniqueneighborvals(uniqueneighborvals==0) = [];
            
            %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
            if isempty(uniqueneighborvals)
                clusteredmetric(vertexneighbors_thiscolor) = vertex;
                %if there is only one previous cluster identifier present, make all the neighbors that value
            elseif length(uniqueneighborvals)==1
                clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
                %if there are multiple cluster identifier values in the neighborhood, merge them into one
            else
                for valuenum = 2:length(uniqueneighborvals)
                    clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
                end
            end
        end
        uniqueclustervals = unique(clusteredmetric);
        uniqueclustervals(uniqueclustervals==0) = [];
        
        for clusternum = uniqueclustervals'
            if nnz(clusteredmetric==clusternum) < sizethresh
                neighborverts = unique(neighbors((clusteredmetric==clusternum),2:7));
                neighborverts(isnan(neighborverts)) = [];
                borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
                borderverts(hemdata(borderverts)<1) = [];
                mode_neighborval = mode(hemdata(borderverts));
                hemdata(clusteredmetric==clusternum) = mode_neighborval;
            end
        end
    end
    
    cleaned([1:nnz(mask)] + (nnz(maskL) * (hemnum-1))) = hemdata(logical(mask));
end

out = zeros(size(groupnetworks));
out(1:nsurfverts) = cleaned;
outfilename = [regularizednetworksfolder '/Consensus_Powercolors_cleaned_' regularizednetworksfile];
cifti_write_wHDR(out,[],outfilename)
    