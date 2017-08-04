hem = 'L';
pialfile = '';
whitefile = '';
vertexnum = 1;
circle_radius = 10;
white_weighting = [.6 .7 .8 .9 .99];
outname_base = '';
load_distances = 1;







if load_distances
    distances = smartload(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);
end

vert_distance = distances(vertexnum);
circle_indices = vert_distance <= circle_radius;

pial = gifti(pialfile);
white = gifti(whitefile);

for weight = white_weighting;
    
    out = pial;
    out.vertices(circle_indices,:) = (pial.vertices(circle_indices,:) * (1-weight)) + (white.vertices(circle_indices,:) * weight);
    
    save(out,[outname_base '_whiteweighted' num2str(weight) '.surf.gii'])
end



