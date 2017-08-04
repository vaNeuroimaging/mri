function make_striped_cifti(subid,hubday,thresh,minclustersizemm,outname,stripewidth)
% make_striped_cifti(inputmaps,minclustersizemm,outname,[stripewidth])
%
% This function takes a multi-column input of cifti maps with categorical,
% non-continous values arranged in clusters on the surface. The function
% collapses across columns, representing regions that have multiple
% identities aross columns with alternating stripes.
%
% Subcortical data is removed (can't easily make stripes in subcortex).
%
% Requires the ft_read/write_cifti_mod scripts to be in your path.
%
% Inputs
%  inputmaps: cifti data (or a string path to cifti data) with categorical
%     values in multiple columns. Must be in fs_LR_32k space with data
%     excluded from the atlas medial wall.
%  minclustersizemm: the minimum size of contiguous single-value clusters
%     in the input that will be represented in the overlap image
%  outname: the name of the output file (with no extension). The output
%     cifti will be in fs_LR_164k space with data excluded from the atlas
%     medial wall. 
%  stripewidth: optional argument. Approximate width of stripes in radians
%     on the spherical expansion of the cortical surface. Defaults to 1/40.
%
% E. Gordon 03/18/15 

%%%%%%%%%%%% begin CL edits %%%%%%%%%%%%

%% define systems
systems = 1:14; 
costs = 90:1:99;
%% load fc matrix, network vector, and hub parcel
load(['/share/cjv2/data/' subid '/import/Day1/Rest/FC_matrix_d1.mat']);
load(['/share/cjv2/data/' subid '/import/Day1/Rest/Networks_d1.mat']);
load(['/share/cjv2/data/' subid '/import/Day1/Rest/HubParcel.mat']);

%% sort networks and parcel X parcel matrix
[~,indx] = sort(net_orig); net = sort(net_orig); matrix = matrix(indx,indx);
matrix_orig_rest(:,:,1) = matrix;

%% loop through costs and construct parcel x parcel x cost matrix
for i = 1:length(costs)
matrix = matrix_orig_rest(:,:,1);
matrix(matrix<prctile(icatb_mat2vec(matrix),costs(i)))=0;
matrix(matrix~=0)=1;
matrix_total(:,:,(i)) = matrix;
end

%% sum accross third dimension
matrix = sum(matrix_total,3);

%% for each node, compute pertcentage connections to each network
for node_i = 1:length(matrix)
for system_ii = 1:length(systems)  
nodes_percent(node_i,system_ii) = sum(matrix(node_i,(net==system_ii)))/sum(matrix(node_i,:));
end 
end

%% find hub row, isolate hub row  
target_row = find(indx==target_row);
hub_percent=nodes_percent((target_row),:);

%% load hub parcel
hub = ft_read_cifti_mod(['/share/cjv2/data/' subid '/import/Day1/targets/Day' num2str(hubday) '/files/' subid '.Target_Day' num2str(hubday) '.dtseries.nii']);
inputmaps = hub; inputmaps.data = zeros(length(inputmaps.data),1);

%% construct inputmap
for ii = 1:length(systems)    

percent = hub_percent(ii);
if percent > thresh
inputmaps.data((find(hub.data~=0)),(ii)) = (ii); 
else
inputmaps.data((find(hub.data~=0)),(ii)) = 0; 
end

end

inputmaps = inputmaps.data;

%%%%%%%%%%%% end of CL edits %%%%%%%%%%%%

outname = num2str(outname);

if ~exist('stripewidth')
    stripewidth = 1/40;
end

spherical_surfaces_32k = {'/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.surf.gii','/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.surf.gii'};
spherical_surfaces_164k = {'/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.L.sphere.164k_fs_LR.surf.gii','/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.R.sphere.164k_fs_LR.surf.gii'};

surfaceareas = ft_read_cifti_mod('/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_cortexonly.dtseries.nii');
inputmaps = inputmaps(1:size(surfaceareas.data,1),:);
neighbors = smartload('/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_surf_neighbors_LR_normalwall.mat');

mask_32k{1} = gifti(['/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/parcellations_VGD11b.L.32k_fs_LR.label.gii']); mask_32k{1} = ~mask_32k{1}.cdata(:,3);
mask_32k{2} = gifti(['/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/parcellations_VGD11b.R.32k_fs_LR.label.gii']); mask_32k{2} = ~mask_32k{2}.cdata(:,3);

mask_164k{1} = gifti(['/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas.LR.164k_fs_LR/parcellations_VGD11b.L.164k_fs_LR.label.gii']); mask_164k{1} = ~mask_164k{1}.cdata(:,3);
mask_164k{2} = gifti(['/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas.LR.164k_fs_LR/parcellations_VGD11b.R.164k_fs_LR.label.gii']); mask_164k{2} = ~mask_164k{2}.cdata(:,3);

cifti_out_templatefile = '/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.LR.164k_fs_LR.normalwall_surfaceonly_template.dtseries.nii';
cifti_out = ft_read_cifti_mod(cifti_out_templatefile);

hems = {'L','R'};

networks = unique(inputmaps); networks(networks==0) = [];

for networknum = 1:length(networks)
    
    networkpresent = any(inputmaps==networks(networknum),2);
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
        
        if sum(surfaceareas.data(clusteredmetric==uniqueclustervals(clusternum))) < minclustersizemm
            
            indicestozero = find(clusteredmetric==uniqueclustervals(clusternum));
            
            for indexnum = indicestozero'
                inputmaps(indexnum,(inputmaps(indexnum,:)==networks(networknum))) = 0;
            end
            
        end
    end
    
end


sorted_inputmaps = sort(inputmaps,2);
[ign, ign2, valuecombinations] = unique(sorted_inputmaps,'rows');

temp = surfaceareas; temp.data = valuecombinations;
ft_write_cifti_mod('Temp',temp);

system(['wb_command -cifti-resample Temp.dtseries.nii COLUMN ' cifti_out_templatefile ' COLUMN BARYCENTRIC ENCLOSING_VOXEL Temp_164.dtseries.nii -surface-largest -left-spheres ' spherical_surfaces_32k{1} ' ' spherical_surfaces_164k{1} ' -right-spheres ' spherical_surfaces_32k{2} ' ' spherical_surfaces_164k{2}]);

valuecombinations_upsampled = ft_read_cifti_mod('Temp_164.dtseries.nii'); valuecombinations_upsampled = valuecombinations_upsampled.data;

surf_withlines = zeros(size(valuecombinations_upsampled));

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    sphere = gifti(['/share/cjv2/software/Software/Gordon_2014_parcellation/parcellation/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.' hem '.sphere.164k_fs_LR.surf.gii']);
    [phi, theta, r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    thetavals = -pi/2 : stripewidth : pi/2;
    surf_withlines_hem = zeros(size(phi));
    
    
    
    for i = 1:length(thetavals)-1
        theseindices = intersect(find(theta>thetavals(i)), find(theta>thetavals(i+1)));
        surf_withlines_hem(theseindices) = i;
    end
    
    surf_withlines([1:nnz(mask_32k{hemnum})] + (hemnum-1)*(nnz(mask_32k{1}))) = surf_withlines_hem(logical(mask_32k{hemnum}));
end

surf_withlines = surf_withlines+1;

finaloutput = zeros(size(valuecombinations_upsampled));
listofcombinations = unique(valuecombinations);
for combinationnum = 1:length(listofcombinations)
    
    combination_withlines = surf_withlines .* (valuecombinations_upsampled==listofcombinations(combinationnum));
    lines_within_parcel = unique(combination_withlines);
    lines_within_parcel(lines_within_parcel==0)=[];
    
    if ~isempty(lines_within_parcel)
        
        values_toassign = sorted_inputmaps(valuecombinations==listofcombinations(combinationnum),:);
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

save('/home/cl968/scripts/inputmaps','inputmaps')

cifti_out.data = finaloutput;
ft_write_cifti_mod(outname,cifti_out)
delete Temp.dtseries.nii Temp_164.dtseries.nii

