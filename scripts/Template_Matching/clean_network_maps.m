function clean_network_maps(networksfile,sizethresh,fillin,subject)

data = ft_read_cifti_mod(networksfile);
neighbors = cifti_neighbors(networksfile);
ncortverts = nnz(data.brainstructure==1) + nnz(data.brainstructure==2);

if exist('subject')
    hems = {'L','R'};
    for hemnum = 1:length(hems)
        surfaceareafiles{hemnum} = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
        surfaceareas{hemnum} = gifti(surfaceareafiles{hemnum});
        surfaceareas{hemnum} = surfaceareas{hemnum}.cdata;
        surfaceareas{hemnum}(data.brainstructure((1:length(surfaceareas{hemnum})) + (length(surfaceareas{hemnum}) * (hemnum-1)))==-1) = [];
    end
    
    surfacearea_voxvol = [surfaceareas{1} ; surfaceareas{2} ; (ones(size(correlmaps_thresh,1) - ncortverts , 1) * voxvol)];
else
    surfacearea_voxvol = ft_read_cifti_mod('/home/data/atlases/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.midthickness.32k_fs_LR_surfaceareas.dtseries.nii');
    surfacearea_voxvol = surfacearea_voxvol.data;
end

temp = zeros(size(data.data,1),1);
IDs = unique(data.data); IDs(IDs<1) = [];
for ID = IDs(:)'
    clustereddata = cifti_cluster_surfacearea_volume(data.data,ID-.5,ID+.5,sizethresh,sizethresh,ncortverts,surfacearea_voxvol,neighbors);
    clustereddata(logical(clustereddata)) = ID;
    temp = temp + clustereddata;

end

data.data = temp;


if exist('fillin') && logical(fillin)
    %-----------------------------
    %Fill in removed spots
    orig = data.data;
    if ~any(orig(ncortverts+1:end))
        data.data(ncortverts+1:end) = 1;
    end
    
    blankinds = find(data.data==0);
    while ~isempty(blankinds)
        temp = data.data;
        for ind = blankinds(:)'
            indneighs = neighbors(ind,2:end); indneighs(isnan(indneighs)) = [];
            neighvals = data.data(indneighs); neighvals(neighvals==0) = [];
            if ~isempty(neighvals)
                temp(ind) = mode(neighvals);
            end
        end
        data.data = temp;
        blankinds = find(data.data==0);
    end
    
    if ~any(orig(ncortverts+1:end))
        data.data(ncortverts+1:end) = 0;
    end
end

dotslocations = strfind(networksfile,'.');
outname = [networksfile(1:(dotslocations(end-1)-1)) '_cleaned'];
ft_write_cifti_mod(outname,data);

