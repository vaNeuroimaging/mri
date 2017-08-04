function make_striped_cifti(inputmaps,minclustersizemm,outname,stripewidth)
% make_striped_cifti(inputmaps,minclustersizemm,outname,[stripewidth])
%
% This function takes a multi-column input of cifti maps with categorical,
% non-continous values arrannged in clusters on the surface. The function
% collapses across columns, representing regions that have multiple
% identities aross columns with alternating stripes.
%
% Requires the ft_read/write_cifti_mod scripts to be in your path.
%
% Inputs
%  inputmaps: cifti data (or a string path to cifti data) with categorical
%     values in multiple columns. Must be in fs_LR_32k space with the atlas
%     medial wall.
%  minclustersizemm: the minimum size of contiguous single-value clusters
%     in the input that will be represented in the overlap image
%  outname: the name of the output file (with no extension). The output
%     cifti will be in fs_LR_164k space with the atlas medial wall.
%  stripewidth: optional argumant. Approximate width of stripes in radians
%     on the spherical expansion of the cortical surface. Defaults to 1/70.
%
% E. Gordon 03/18/15

if ~exist('stripewidth')
    stripewidth = 1/70;
end

surfaceareas = ft_read_cifti_mod('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii'); %surfaceareas = surfaceareas.data;

if ischar(inputmaps)
    inputmaps = ft_read_cifti_mod(inputmaps); inputmaps = inputmaps.data(1:size(surfaceareas.data,1),:);
end

neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');

mask_32k{1} = gifti(['/data/cn4/evan/fsaverage_LR32k/parcellations_VGD11b.L.32k_fs_LR.label.gii']); mask_32k{1} = ~mask_32k{1}.cdata(:,3);
mask_32k{2} = gifti(['/data/cn4/evan/fsaverage_LR32k/parcellations_VGD11b.R.32k_fs_LR.label.gii']); mask_32k{2} = ~mask_32k{2}.cdata(:,3);

mask_164k{1} = gifti(['/data/cn4/evan/Conte_164k/parcellations_VGD11b.L.164k_fs_LR.label.gii']); mask_164k{1} = ~mask_164k{1}.cdata(:,3);
mask_164k{2} = gifti(['/data/cn4/evan/Conte_164k/parcellations_VGD11b.R.164k_fs_LR.label.gii']); mask_164k{2} = ~mask_164k{2}.cdata(:,3);

cifti_out_templatefile = '/data/cn4/evan/Conte_164k/Conte69.LR.164k_fs_LR.normalwall_surfaceonly_template.dtseries.nii';
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
[combinations, ign, valuecombinations] = unique(sorted_inputmaps,'rows');

temp = surfaceareas; temp.data = valuecombinations;
ft_write_cifti_mod('Temp',temp);

system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -cifti-resample Temp.dtseries.nii COLUMN ' cifti_out_templatefile ' COLUMN BARYCENTRIC ENCLOSING_VOXEL Temp_164.dtseries.nii -surface-largest -left-spheres /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.L.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.L.sphere.164k_fs_LR.surf.gii -right-spheres /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.R.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.R.sphere.164k_fs_LR.surf.gii'])

valuecombinations_upsampled = ft_read_cifti_mod('Temp_164.dtseries.nii'); valuecombinations_upsampled = valuecombinations_upsampled.data;

surf_withlines = zeros(size(valuecombinations_upsampled));

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    sphere = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii']);
    [phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    thetavals = -pi/2 : stripewidth : pi/2;
    surf_withlines_hem = zeros(size(phi));
    
    for i = 1:length(thetavals)-1
        theseindices = intersect(find(theta>thetavals(i)), find(theta>thetavals(i+1)));
        surf_withlines_hem(theseindices) = i;
    end
    
    surf_withlines([1:nnz(mask_164k{hemnum})] + (hemnum-1)*(nnz(mask_164k{1}))) = surf_withlines_hem(logical(mask_164k{hemnum}));
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

cifti_out.data = finaloutput;
ft_write_cifti_mod(outname,cifti_out)
delete Temp.dtseries.nii Temp_164.dtseries.nii

% 
% for hemnum = 1:length(hems)
%     hem = hems{hemnum};
%     
%     sorted_inputmaps_32k = zeros(size(mask_32k{hemnum},1),size(sorted_inputmaps,2));
%     sorted_inputmaps_32k(logical(mask_32k{hemnum}),:) = sorted_inputmaps(([1:nnz(mask_32k{hemnum})] + (hemnum-1)*(nnz(mask_32k{1}))),:);
%     
%     valuecombinations_32k = zeros(size(mask_32k{hemnum}));
%     valuecombinations_32k(logical(mask_32k{hemnum})) = valuecombinations(([1:nnz(mask_32k{hemnum})] + (hemnum-1)*(nnz(mask_32k{1}))));
%     
%     save(gifti(single(valuecombinations_32k)),'Temp.func.gii');
%     evalc(['!/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Temp.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Temp_164.func.gii -largest']);
%     valuecombinations_upsampled = gifti('Temp_164.func.gii'); valuecombinations_upsampled = valuecombinations_upsampled.cdata;
%     
%     sphere = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii']);
%     [phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
%     thetavals = -pi/2 : stripewidth : pi/2;
%     surf_withlines = zeros(size(phi));
%     
%     for i = 1:length(thetavals)-1
%         theseindices = intersect(find(theta>thetavals(i)), find(theta>thetavals(i+1)));
%         surf_withlines(theseindices) = i;
%     end
%     
%     surf_withlines = surf_withlines+1;
%     
%     finaloutput = zeros(size(valuecombinations_upsampled));
%     listofcombinations = unique(valuecombinations_32k);
%     for combinationnum = 1:length(listofcombinations)
%         
%         combination_withlines = surf_withlines .* (valuecombinations_upsampled==listofcombinations(combinationnum));
%         lines_within_parcel = unique(combination_withlines);
%         lines_within_parcel(lines_within_parcel==0)=[];
%         
%         if ~isempty(lines_within_parcel)
%             
%             values_toassign = sorted_inputmaps_32k(valuecombinations_32k==listofcombinations(combinationnum),:);
%             values_toassign = values_toassign(1,:);
%             values_toassign(values_toassign==0) = [];
%             
%             
%             if length(values_toassign) <= length(lines_within_parcel)
%                 
%                 for comnum = 1:length(values_toassign)
%                     
%                     lineval_indices = [comnum : length(values_toassign) : length(lines_within_parcel)];
%                     linevals_for_this_community = lines_within_parcel(lineval_indices);
%                     for lineval = linevals_for_this_community'
%                         finaloutput(combination_withlines==lineval) = values_toassign(comnum);
%                     end
%                     
%                 end
%             end
%         end
%     end
%     
%     cifti_out.data([1:nnz(mask_164k{hemnum})] + (hemnum-1)*(nnz(mask_164k{1}))) = finaloutput(logical(mask_164k{hemnum}));
%     
% end
%
% ft_write_cifti_mod(outname,cifti_out)
% delete Temp.func.gii Temp_164.func.gii
