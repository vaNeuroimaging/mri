function consensus_maker_knowncolors_yeo(regularized_ciftifile,groupnetworksfile,mincol,minsize,recolor)

if ~exist('mincol','var') || isempty(mincol)
    mincol = 1;
end

if ~exist('minsize','var') || isempty(minsize)
    minsize = 0;
end

if ~exist('groupnetworksfile','var') || isempty(groupnetworksfile)
    groupnetworksfile = '/home/data/evan/Temp/Yeo_7networks_powercolors.dtseries.nii';
end

% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
all_color_values = [1:17 2.5 4.5 11.5 18:100];

cifti_data = ft_read_cifti_mod(regularized_ciftifile); assigns = cifti_data.data;
assigns(assigns<0) = 0;
assigns(isnan(assigns)) = 0;




groupfile = ft_read_cifti_mod(groupnetworksfile);
groupdata = groupfile.data;
groupdata = groupdata(1:59412); groupdata(59413:size(assigns,1)) = 0;
%potential_colors = unique(groupdata); potential_colors(potential_colors<=0) = [];
potential_colors = [1 2 10 9 5 3 11 16 15 7 8 12 14 13]; 
newcolors = setdiff(all_color_values,potential_colors);
%colors(colors==0) = [];


all_recolored = zeros(size(assigns));
for c = 1:size(all_recolored,2)
    col_consensusmap = assigns(:,c);
    
    
    unassigned = find(col_consensusmap<1);
    for unassignedindex = unassigned'
        thisassignments = assigns(unassignedindex,mincol:end);
        thisassignments(thisassignments<1) = [];
        if ~isempty(thisassignments)
            col_consensusmap(unassignedindex) = thisassignments(1);
        end
    end
    networks = unique(col_consensusmap); networks(networks<=0) = [];
    new_networks = networks;
    assigning_networks = networks;
    
    col_out = zeros(size(col_consensusmap));
    for i = 1:length(potential_colors)
        if ~isempty(assigning_networks)
            groupnetwork_comp = groupdata==potential_colors(i);
            D = zeros(length(assigning_networks),1);
            P = zeros(length(assigning_networks),1);
            for j = 1:length(assigning_networks)
                
                network_comp = col_consensusmap==assigning_networks(j);
                P(j) = nnz(groupnetwork_comp & network_comp);
                D(j) = P(j)/nnz(groupnetwork_comp | network_comp);
            end
                [maxval, maxind(i)] = max(D);
                %[maxval maxind_P(i)] = max(P);
            if maxval > .1
                col_out(col_consensusmap==assigning_networks(maxind(i))) = potential_colors(i);
                new_networks(new_networks==assigning_networks(maxind(i))) = [];
                assigning_networks(assigning_networks==assigning_networks(maxind(i))) = [];
            end
        end
        
    end
    clear maxind D P
    for j = 1:length(new_networks)
        col_out(col_consensusmap==new_networks(j)) = newcolors(j);
    end
    all_recolored(:,c) = col_out;
end

out = all_recolored(:,mincol);

uniquevals = unique(out); uniquevals(uniquevals<1) = [];
colors_tofix = setdiff(uniquevals,potential_colors);
verts_tofix = [];
for colornum = 1:length(colors_tofix)
    verts_thiscolor = find(out==colors_tofix(colornum));
    verts_tofix = [verts_tofix ; verts_thiscolor(:)];
end
for vertnum = verts_tofix'
    for col = (mincol+1):size(all_recolored,2)
        if any(all_recolored(vertnum,col)==potential_colors)
            out(vertnum) = all_recolored(vertnum,col);
            break
        end
    end
end




if logical(minsize)


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

end

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
outname = [basename '_recolored_yeo'];
ft_write_cifti_mod(outname,cifti_data);



        
    
