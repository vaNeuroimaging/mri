function rotate_values_on_sphere(metricname,rotation,hem)

metric = gifti(metricname);
metric = metric.cdata;


clear indpos
surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);

[phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));


phi_new = phi+pi;
theta_new = theta+pi/2;

phi_new = mod(phi_new+rotation,2*pi)-pi;



for n = 1:length(phi_new)
    
    test = repmat([phi_new(n) theta(n)],[32492 1]);
    diff_test = sum(abs([phi theta] - test),2);
    [val ind(n)] = min(diff_test);
end

outputmetric = metric(ind);

save(gifti(single(outputmetric)),[metricname(1:end-9) '_rotated' num2str(rotation) '.func.gii'])






