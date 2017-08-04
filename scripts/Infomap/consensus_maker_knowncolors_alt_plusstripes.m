function consensus_maker_knowncolors_alt_plusstripes(regularized_ciftifile,manualset,groupnetworksfile,mincol,minsize,orig_parcelsfile,recolor)
%consensus_maker_knowncolors_alt_plusstripes(regularized_ciftifile,manualset,groupnetworksfile,mincol,minsize,orig_parcelsfile,recolor)

if ~exist('mincol','var') || isempty(mincol)
    mincol = 1;
end

if ~exist('minsize','var') || isempty(minsize)
    minsize = 0;
end

if ~exist('groupnetworksfile','var') || isempty(groupnetworksfile)
    groupnetworksfile = '/home/data/atlases/120_LR_minsize400_recolored_manualconsensus4.dtseries.nii';
end

% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
all_color_values = [1:100];

cifti_data = ft_read_cifti_mod(regularized_ciftifile); assigns = cifti_data.data;
assigns(assigns<0) = 0;
assigns(isnan(assigns)) = 0;




groupfile = ft_read_cifti_mod(groupnetworksfile);
groupdata = groupfile.data;
if size(groupdata,1) ~= size(assigns,1)
    groupdata = groupdata(1:59412);
    assigns = assigns(1:59412,:);
else
    groupdata(59413:size(assigns,1)) = 0;% groupdata = groupdata(1:min(59412,size(groupdata,1))); 
end
%potential_colors = unique(groupdata); potential_colors(potential_colors<=0) = [];
potential_colors = [1 2 10 9 3 5 6 11 16 15 7 8 17 12 4 14 13]; 
newcolors = setdiff(all_color_values,potential_colors);
%colors(colors==0) = [];

unassigned_networks = cell(1,size(assigns,2));

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
    
    if exist('manualset')
        for i = 1:size(manualset,1)
            col_out(col_consensusmap==manualset(i,1)) = manualset(i,2);
            new_networks(new_networks==manualset(i,1)) = [];
            if all(manualset(i,2)~=potential_colors); new_networks = [new_networks; manualset(i,2)]; end
            assigning_networks(assigning_networks==manualset(i,1)) = [];
        end
    end
    
    
    for i = 1:length(potential_colors)
        
        if exist('manualset') && any(manualset(:,2)==potential_colors(i))
            
%             manualassignedinds = find(manualset(:,2)==potential_colors(i));
%             for manualassignedind = manualassignedinds(:)'
%                 col_out(col_consensusmap==manualset(manualassignedind,1)) = potential_colors(i);
%                 new_networks(new_networks==manualset(manualassignedind,1)) = [];
%                 assigning_networks(assigning_networks==manualset(manualassignedind,1)) = [];
%             end
            
                
        else

        
        if ~isempty(assigning_networks)
            groupnetwork_comp = groupdata==potential_colors(i);
            D = zeros(length(assigning_networks),1);
            P = zeros(length(assigning_networks),1);
            for j = 1:length(assigning_networks)
                
                %if c>1 && any(unassigned_networks{c-1}==assigning_networks(j)) && any(any(all_recolored(:,1:(c-1))==potential_colors(i)))
                %    D(j) = 0;
                %else
                    network_comp = col_consensusmap==assigning_networks(j);
                    P(j) = nnz(groupnetwork_comp & network_comp);
                    D(j) = P(j)/nnz(groupnetwork_comp | network_comp);
                %end
            end
                [maxval, maxind(i)] = max(D);
                if potential_colors(i)==17 || potential_colors(i)==14;
                1;
                end
            if maxval > .1

                col_out(col_consensusmap==assigning_networks(maxind(i))) = potential_colors(i);
                new_networks(new_networks==assigning_networks(maxind(i))) = [];
                assigning_networks(assigning_networks==assigning_networks(maxind(i))) = [];
            end
        end
        end
        
    end
    clear maxind D P
    for j = 1:length(new_networks)
        col_out(col_consensusmap==new_networks(j)) = newcolors(j);
    end
    all_recolored(:,c) = col_out;
    unassigned_networks{c} = assigning_networks;
end


all_recolored(assigns<=0) = 0;
cifti_data.data(1:size(all_recolored,1),:) = all_recolored;
if ~exist('cifti_data.mapname')
    for col = 1:size(all_recolored,2)
    cifti_data.mapname{col} = ['Column number ' num2str(col)];
    end
    cifti_data.dimord = 'scalar_pos';
end

dotsloc = strfind(regularized_ciftifile,'.');
basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_allcolumns_recoloredv4'];
ft_write_cifti_mod(outname,cifti_data);
set_cifti_powercolors([outname '.dscalar.nii'])




out = all_recolored(:,mincol);

uniquevals = unique(out); %uniquevals(uniquevals<1) = [];
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


temp_out = out;

if logical(minsize)


% Clean up tiny pieces


neighbors = smartload('/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_neighbors_normalwall.mat');

allcolors = unique(out); allcolors(allcolors<=0) = [];


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




cifti_data.data = cifti_data.data(:,1);
cifti_data.data(1:size(out,1),:) = out;
if ~exist('cifti_data.mapname')
    cifti_data.mapname = {'Column number ' num2str(mincol)};
    cifti_data.dimord = 'scalar_pos';
else
    cifti_data.mapname = cifti_data.mapname(mincol);
end

dotsloc = strfind(regularized_ciftifile,'.');
basename = regularized_ciftifile(1:(dotsloc(end-1)-1));
outname = [basename '_recoloredv4'];
ft_write_cifti_mod(outname,cifti_data);

set_cifti_powercolors([outname '.dscalar.nii'])

if exist('orig_parcelsfile') && ~isempty(orig_parcelsfile)
    parcels = ft_read_cifti_mod(orig_parcelsfile);
    parcels = parcels.data;
    parcels((length(out)+1):end) = [];
    IDs = unique(parcels); IDs(IDs<1) = [];
    outtext = zeros(length(IDs),1);
    for IDnum = 1:length(IDs)
        outtext(IDnum) = mode(out(parcels==IDs(IDnum)));
        
    end
    dlmwrite([outname '.txt'],outtext,'delimiter',' ')
end

      
%Stripes

all_recolored = all_recolored(1:59412,:);
allcolors = unique(all_recolored); allcolors(allcolors==0) = [];
unknown_colors = setdiff(allcolors,potential_colors);
for color = unknown_colors(:)'
    all_recolored(all_recolored==color) = 0;
end

to_be_striped = out;
change = logical(diff(all_recolored,1,2));
change = change .* (all_recolored(:,1:end-1)>0) .* (all_recolored(:,2:end)>0);
for col = 1:size(change,2)
    colvals = unique(all_recolored(:,col+1)); colvals(colvals==0) = [];
    for val_totest = colvals(:)'
        if nnz((all_recolored(:,col)==val_totest) & (all_recolored(:,col+1)==val_totest)) / nnz((all_recolored(:,col)==val_totest) | (all_recolored(:,col+1)==val_totest)) > .5
            thiscolor_changed_inds = logical(change(:,col) .* (all_recolored(:,col)==val_totest));
            to_be_striped(thiscolor_changed_inds,2+((col-1)*2)) = all_recolored(thiscolor_changed_inds,col+1);
            to_be_striped(thiscolor_changed_inds,3+((col-1)*2)) = val_totest;
        end
    end
end

        
to_be_striped_final = zeros(size(to_be_striped));
for vert = 1:size(to_be_striped_final,1)
    uniquevals = unique([out(vert) to_be_striped(vert,:)]); uniquevals(uniquevals==0) = [];
    to_be_striped_final(vert,1:length(uniquevals)) = uniquevals;
end
make_striped_cifti(to_be_striped_final,0,[outname '_striped_164'],1/40)
set_cifti_powercolors([outname '_striped_164.dtseries.nii'])    
