function ind = rotate_values_on_sphere_xyz_func(num,hem)   

surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
spherecoords = sphere.vertices;

for r = 1:num
    xrot = rand * 2*pi;
    yrot = rand * 2*pi;
    zrot = rand * 2*pi;
    
    rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
    rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
    rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
    
    indexcoords = spherecoords';
    xrotcoords = rotmat_x * indexcoords;
    xyrotcoords = rotmat_y * xrotcoords;
    xyzrotcoords = rotmat_z * xyrotcoords;
    
    
    
    for n = 1:size(xyzrotcoords,2)
        testcoord = xyzrotcoords(:,n)';
        test = repmat(testcoord,[32492 1]);
        diff_test = sum(sqrt((spherecoords - test).^2),2);
        [val ind(n,r)] = min(diff_test);
    end
end