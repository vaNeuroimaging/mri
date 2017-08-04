function neighbors = cifti_neighbors_nostructdivisions(ciftifile)
%neighbors = cifti_neighbors(ciftifile)

if isstruct(ciftifile)
    cifti = ciftifile; clear ciftifile;
else
    cifti = ft_read_cifti_mod(ciftifile); 
end
    cifti.data = [];

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[surfneighbors(:,1) surfneighbors(:,2) surfneighbors(:,3) surfneighbors(:,4)...
    surfneighbors(:,5) surfneighbors(:,6) surfneighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
surfneighbors = surfneighbors+1;

nsurfverts = size(surfneighbors,1);

for hem = 1:2
    
    hem_surfneighbors{hem} = surfneighbors;
    verts_notindata = find(cifti.brainstructure((1:nsurfverts)+(nsurfverts*(hem-1)))==-1);
    hem_surfneighbors{hem}(verts_notindata,:) = [];
    for vert_notindata = verts_notindata(:)'
        hem_surfneighbors{hem}(hem_surfneighbors{hem}==vert_notindata) = NaN;
    end
    
    temp = hem_surfneighbors{hem};
    
    verts_indata = find(cifti.brainstructure((1:nsurfverts)+(nsurfverts*(hem-1)))==hem);
    for vertnum_indata = 1:length(verts_indata)
        vert_indata = verts_indata(vertnum_indata);
        temp(hem_surfneighbors{hem}==vert_indata) = vertnum_indata + (nnz(cifti.brainstructure==1)*(hem-1));
    end
    hem_surfneighbors{hem} = temp;
        
    hem_surfneighbors{hem}(:,end+1:27) = NaN;
    
end

neighbors = [hem_surfneighbors{1} ; hem_surfneighbors{2}];


if isfield(cifti,'transform')
    
    dims = [cifti.transform(1,1) cifti.transform(2,2) cifti.transform(3,3)];
    toofar_tobeneighbors_dist = min([sqrt(dims(1)^2 + dims(2)^2) sqrt(dims(1)^2 + dims(3)^2) sqrt(dims(2)^2 + dims(3)^2)]);
    volpos = cifti.pos(cifti.brainstructure>2,:);
    volneighbors = ones(size(volpos,1),27) * NaN;
    for vox = 1:size(volpos,1)
        distances_fromthisvox = pdist2(volpos(vox,:),volpos);
        thisvox_neighbors = find(distances_fromthisvox < toofar_tobeneighbors_dist);
        thisvox_neighbors(thisvox_neighbors==vox) = [];
        thisvox_neighbors = [vox; thisvox_neighbors(:)];
        volneighbors(vox,1:length(thisvox_neighbors)) = thisvox_neighbors + nnz(cifti.brainstructure==1) + nnz(cifti.brainstructure==2);
    end
    neighbors = [neighbors; volneighbors];
else
    neighbors = neighbors(:,1:7);
end

