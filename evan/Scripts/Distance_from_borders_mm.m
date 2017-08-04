function distborders = Distance_from_borders_mm(threshedges1,threshedges2,hem,outputname)
%distborders = Distance_from_borders_mm(threshedges1,threshedges2,hem,[outputname])

if ischar(threshedges1)
    threshedges1 = gifti(threshedges1); threshedges1 = threshedges1.cdata;
end

if ischar(threshedges2)
    threshedges1 = gifti(threshedges2); threshedges2 = threshedges2.cdata;
end

load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])

distborders = ones(size(threshedges1)) * -1;

for i = find(threshedges1)'
    
        dist = min(geo_distances(i,logical(threshedges2)));
        
        distborders(i) = dist;
            
end

if exist('outputname')
    save(gifti(single(distborders)),outputname);
end