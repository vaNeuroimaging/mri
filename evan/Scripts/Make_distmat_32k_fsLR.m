function Make_distmat_32k_fsLR(ciftifile,outputfile,volume_space)
%Make_distmat_32k_fsLR(ciftifile,outputfile,[volume_space])
%
% Given a cifti in 32k_fsLR space (with any subcortical components), this
% function will create a point-to-point distance matrix and save it to a
% specified output file. This matrix is useful for many functions,
% particularly infomap.
%
% The distance matrix will use geodesic distances for the distance from each
% cortical point to each other cortical point in the same cortical
% hemisphere. Euclidean distance will be used for distances between
% hemispheres, as well as for distances from cortical to volumetric (i.e.
% subcortical) points.
%
% Subcortical points are assumed to be in 711-2b space. If these points are in
% MNI space, the variable 'volume_space' should be specified as 'MNI'.
%
% Inputs:
% ciftifile - the full path to a cifti file in the same space as the
%  desired distance matrix
% outputfile - the full path of the desired output file, including .mat
%  extension
% volume_space - an optional argument specifying the space of the
%  volumetric subcortical elements of the cifti. Specify '711-2b' or 'MNI'.
%  Defaults to '711-2b'.
%
% EMG 08/20/15

surfcoordsL = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.coord.gii'); 
surfcoordsR = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.midthickness.32k_fs_LR.coord.gii'); 
surfcoordsL = surfcoordsL.vertices;
surfcoordsR = surfcoordsR.vertices;

cifti = ft_read_cifti_mod(ciftifile);
maskL = cifti.brainstructure(1:length(surfcoordsL)) > 0;
maskR = cifti.brainstructure((length(surfcoordsL)+1):(length(surfcoordsL)+length(surfcoordsR))) > 0;

surfcoordsL = surfcoordsL(logical(maskL),:);
surfcoordsR = surfcoordsR(logical(maskR),:);

subcort_coords = cifti.pos((length(maskL)+length(maskR))+1:end,:);

if ~exist('volume_space','var') || strcmp(volume_space,'711-2b')
    
    T4_711_to_MNI_mat = [0.953903 -0.003872 -0.021291; -0.010402  0.950932 -0.016412; 0.020852  0.055067  0.946383];
    for i = 1:size(subcort_coords,1)
        subcort_coords(i,:) = subcort_coords(i,:) * T4_711_to_MNI_mat;
    end
    
elseif ~strcmp(volume_space,'MNI')
    error('volume_space must be ''711-2b'' or ''MNI''')
end

all_MNI_coords = [surfcoordsL ; surfcoordsR ; subcort_coords];
distmat_use = single(squareform(pdist(all_MNI_coords)));

distancesL = smartload('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Surface_distances_L.mat'); 
distancesL = single(distancesL(logical(maskL),logical(maskL)));
distancesR = smartload('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Surface_distances_R.mat'); 
distancesR = single(distancesR(logical(maskR),logical(maskR)));

distmat_use(1:length(surfcoordsL),1:length(surfcoordsL)) = distancesL;
distmat_use((length(surfcoordsL)+1) : (length(surfcoordsL) + length(surfcoordsR)),(length(surfcoordsL)+1) : (length(surfcoordsL) + length(surfcoordsR))) = distancesR;


save(outputfile,'distmat_use','-v7.3')
distmat_use = uint8(distmat_use);
[folder filename ext] = fileparts(outputfile);
save([folder '/' filename '_uint8.' ext],'distmat_use','-v7.3')


