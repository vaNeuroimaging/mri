function ind = rotate_values_on_sphere_xyz_cifti(num,cifti_template)

surfcoordsL = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.coord.gii'); 
surfcoordsR = gifti('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.coord.gii'); 
surfcoordsL = surfcoordsL.vertices;
surfcoordsR = surfcoordsR.vertices;

if ~exist('cifti_template')
    cifti_template = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_normalwall_vertexindices.dtseries.nii';
end
maskL = logical(cifti_template.brainstructure(1:size(surfcoordsL,1)) > 0);
maskR = logical(cifti_template.brainstructure((size(surfcoordsL,1)+1) : (size(surfcoordsL,1)+size(surfcoordsR,1))) > 0);

surfcoordsL = surfcoordsL(maskL,:);
surfcoordsR = surfcoordsR(maskR,:);


indL = zeros(size(surfcoordsL,1),num);
for r = 1:num
    xrot = rand * 2*pi;
    yrot = rand * 2*pi;
    zrot = rand * 2*pi;
    
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    indexcoords = surfcoordsL';
    xrotcoords = rotmat_x * indexcoords;
    xyrotcoords = rotmat_y * xrotcoords;
    xyzrotcoords = rotmat_z * xyrotcoords;
    
    
    
    for n = 1:size(xyzrotcoords,2)
        testcoord = xyzrotcoords(:,n)';
        test = repmat(testcoord,[32492 1]);
        diff_test = sum(sqrt((surfcoordsL - test).^2),2);
        [~,indL(n,r)] = min(diff_test);
    end
end

indR = zeros(size(surfcoordsR,1),num);
for r = 1:num
    xrot = rand * 2*pi;
    yrot = rand * 2*pi;
    zrot = rand * 2*pi;
    
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    indexcoords = surfcoordsR';
    xrotcoords = rotmat_x * indexcoords;
    xyrotcoords = rotmat_y * xrotcoords;
    xyzrotcoords = rotmat_z * xyrotcoords;
    
    
    
    for n = 1:size(xyzrotcoords,2)
        testcoord = xyzrotcoords(:,n)';
        test = repmat(testcoord,[32492 1]);
        diff_test = sum(sqrt((surfcoordsR - test).^2),2);
        [~,indR(n,r)] = min(diff_test);
    end
end

ind = [indL ; (indR + size(indL,1))];