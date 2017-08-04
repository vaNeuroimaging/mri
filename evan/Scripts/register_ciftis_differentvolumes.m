function [outcifti1, outcifti2] = register_ciftis_differentvolumes(incifti1, incifti2)
%[outcifti1, outcifti2] = register_ciftis_differentvolumes(incifti1, incifti2)


%Load data

if ischar(incifti1)
    incifti1 = ft_read_cifti_mod(incifti1);
end

if ischar(incifti2)
    incifti2 = ft_read_cifti_mod(incifti2);
end

if incifti1.transform ~= incifti2.transform
    error('Cifti volume elements not in the same space!')
end




%Set up variables
cifti1_vol_coords = incifti1.pos(incifti1.brainstructure > 2,:);
cifti2_vol_coords = incifti2.pos(incifti2.brainstructure > 2,:);

cifti1_volinds = false(size(cifti1_vol_coords,1),1);
cifti2_volinds = false(size(cifti2_vol_coords,1),1);




%Compare voxel coords
for i = 1:length(cifti1_volinds)
    matchind = ~logical(pdist2(cifti1_vol_coords(i,:),cifti2_vol_coords));
    if any(matchind)
        cifti1_volinds(i) = 1;
        cifti2_volinds(matchind) = 1;
    end
end




%Restrict output data to overlapping voxels: Cifti #1
cifti1_allinds = [true(nnz(incifti1.brainstructure <= 2),1) ; cifti1_volinds];

cifti1_allinds_nomed = cifti1_allinds;
cifti1_allinds_nomed(incifti1.brainstructure<0) = [];

outcifti1 = incifti1; clear incifti1
outcifti1.brainstructure = outcifti1.brainstructure(cifti1_allinds);
outcifti1.pos = outcifti1.pos(cifti1_allinds,:);

if strcmp(outcifti1.dimord,'pos_pos')
    outcifti1.data = outcifti1.data(cifti1_allinds_nomed,cifti1_allinds_nomed);
else
    outcifti1.data = outcifti1.data(cifti1_allinds_nomed,:);
end




%Restrict output data to overlapping voxels: Cifti #1
cifti2_allinds = [true(nnz(incifti2.brainstructure <= 2),1) ; cifti2_volinds];

cifti2_allinds_nomed = cifti2_allinds;
cifti2_allinds_nomed(incifti2.brainstructure<0) = [];

outcifti2 = incifti2; clear incifti2
outcifti2.brainstructure = outcifti2.brainstructure(cifti2_allinds);
outcifti2.pos = outcifti2.pos(cifti2_allinds,:);

% %make brain structure IDs match across datasets
% outcifti2.brainstructure = outcifti1.brainstructure; 
% outcifti2.brainstructurelabel = outcifti1.brainstructurelabel;

if strcmp(outcifti2.dimord,'pos_pos')
    outcifti2.data = outcifti2.data(cifti2_allinds_nomed,cifti2_allinds_nomed);
else
    outcifti2.data = outcifti2.data(cifti2_allinds_nomed,:);
end