function consensus_maker(regularized_ciftifile,mincol,recolor)
%consensus_maker(regularized_ciftifile,mincol,recolor)

if ~exist('mincol','var') || isempty(mincol)
    mincol = 1;
end



% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
all_color_values = [1:17 2.5 4.5 11.5 18:100];
minsize = 15;

cifti_data = ft_read_cifti_mod(regularized_ciftifile); assigns = cifti_data.data;
assigns(assigns<0) = 0;
assigns(isnan(assigns)) = 0;

consensusmap = assigns(:,mincol);


unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = assigns(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end







groupnetworksfile = '/home/data/atlases/120_LR_minsize400_recolored_manualconsensus_LR.dtseries.nii';
groupfile = ft_read_cifti_mod(groupnetworksfile);
groupdata = groupfile.data;
groupdata = groupdata(1:59412); groupdata(59413:length(consensusmap)) = 0;
%potential_colors = unique(groupdata); potential_colors(potential_colors<=0) = [];
potential_colors = [1 2 10 9 5 3 11 16 15 7 8 12 14 13]; 
newcolors = setdiff(all_color_values,potential_colors);

%colors(colors==0) = [];
networks = unique(consensusmap); networks(networks<=0) = [];
new_networks = networks;
assigning_networks = networks;

out = zeros(size(consensusmap));
for i = 1:length(potential_colors)
    
    groupnetwork_comp = groupdata==potential_colors(i);
    
    for j = 1:length(assigning_networks)
        
        network_comp = consensusmap==assigning_networks(j);
        P(j) = nnz(groupnetwork_comp & network_comp);
        D(j) = 2*P(j)/nnz(groupnetwork_comp | network_comp);
    end
    
    [maxval maxind(i)] = max(D);
    %[maxval maxind_P(i)] = max(P);
    if maxval > .05
    out(consensusmap==assigning_networks(maxind(i))) = potential_colors(i);
    new_networks(new_networks==assigning_networks(maxind(i))) = [];
    assigning_networks(assigning_networks==assigning_networks(maxind(i))) = [];
    end
    
end

for i = 1:length(new_networks)
    out(consensusmap==new_networks(i)) = newcolors(i);
end



% Clean up tiny pieces


neighbors = smartload('/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_neighbors_normalwall.mat');

allcolors = unique(out); allcolors(allcolors<=0) = [];

temp_out = out;

for color = allcolors(:)'
    clusteredmetric = zeros(size(temp_out));
    thiscolorverts = find(temp_out==color);
    for vertex = thiscolorverts'
        
        %find the neighbors of this vertex
        vertexneighbors = neighbors(vertex,:);
        vertexneighbors(isnan(vertexneighbors)) = [];
        
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
        if any(find(clusteredmetric==clusternum)==3137)
            1;
        end
        if nnz(clusteredmetric==clusternum) < minsize
            neighborverts = unique(neighbors((clusteredmetric==clusternum),2:end));
            neighborverts(isnan(neighborverts)) = [];
            borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
            borderverts(temp_out(borderverts)<1) = [];
            mode_neighborval = mode(temp_out(borderverts));
            temp_out(clusteredmetric==clusternum) = mode_neighborval;
        end
    end
end


out = temp_out;


if exist('recolor','var')
    for i = 1:size(recolor,1)
        out(temp_out==recolor(i,1)) = recolor(i,2);
    end
end




cifti_data.data = out;
cifti_data.mapname = cifti_data.mapname(mincol);
%cifti_data = rmfield(cifti_data,'mapname');
%cifti_data.dimord = 'pos_time';

dotsloc = strfind(regularized_ciftifile,'.');
basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_recolored'];
ft_write_cifti_mod(outname,cifti_data);



        
    
