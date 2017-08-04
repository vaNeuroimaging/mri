inds = [5 8 11:13 21 24 26 35 38 40:42 43 44 53 54 56 65 66 78 81 83:88 90 91 95 97 100:102 106 111 113 117 122 126:128 130 131 148 152 153 161 163 166 168 178 179 187 192 193 197 200 210:216 218 219 224 261:266 271 278 282 283 286 291];

mapname = 'Cluster_probability_maps_sorted_10mm_40sub';
sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub';
interestingfeature = 'notintemplate_plusextended';

networksfound = {'Default','Vis','FP','DA','VA','Sal','CO','MH','MM','Aud','PERN','ParOcc'};


removegroupcluster_inds = [43 44 95 97 148 278;3 3 7 7 9 16];
dilate = 14;
%neighbors = cifti_neighbors([mapname '.dscalar.nii']);

groupsystems = ft_read_cifti_mod('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
groupsystems = groupsystems.data;

baddata = ft_read_cifti_mod('/data/cn4/evan/ROIs/Baddata_bigcluster_LR_fillholes.dtseries.nii');
baddata = logical(baddata.data);

%distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat');
radius = 5;



maps = ft_read_cifti_mod([mapname '.dscalar.nii']);
sizes = ft_read_cifti_mod([sizesname '.dscalar.nii']);
maps.data = maps.data(:,inds);
sizes.data = sizes.data(:,inds);
maps.mapname = maps.mapname(inds);
sizes.mapname = sizes.mapname(inds);
maps2 = zeros(size(maps.data,1),1);
sizes2 = zeros(size(maps.data,1),1);
col = 1;
tokens = tokenize(maps.mapname{1},' ');
prev = tokens{1};
for i = 1:size(maps.data,2)
    
tokens = tokenize(maps.mapname{i},' ');
thismap = maps.data(:,i) * str2num(tokens{6}) / str2num(tokens{4});
thissizes = sizes.data(:,i);

removenum = find(removegroupcluster_inds(1,:)==inds(i));
if removenum == 5
    1;
end
if ~isempty(removenum)
    ID = removegroupcluster_inds(2,removenum);
    outputclusters = cifti_cluster(groupsystems,ID-.5,ID+.5,0);
    dice = zeros(size(outputclusters,2),1);
    for i = 1:size(outputclusters,2)
        dice(i) = nnz(logical(outputclusters(:,i)) & (thismap>.05)) ./ nnz(logical(outputclusters(:,i)) | (thismap>.05));
    end
    [~,maxdiceind] = max(dice);
    cluster_toremove = outputclusters(:,maxdiceind);
    
    removeinds = find(cluster_toremove);
    dilate_distances = (distances(:,removeinds) < dilate);
    removeinds_dilated = logical(any(dilate_distances,2));
    
    dilate_distances = (distances(:,removeinds) < (dilate+5));
    removeinds_dilated_bigger = logical(any(dilate_distances,2));

    
    thismap(removeinds_dilated) = 0;
    thismap(baddata) = 0;
    
    tempmap = thismap; tempmap(removeinds_dilated_bigger) = 0;
    
    [ign,maxi] = max(tempmap);
    
    thismap(distances(:,maxi) > 50) = 0;
    
    temp = distances(:,maxi) < radius;
    tempneighs = unique(neighbors(logical(temp),2:end));
    tempneighs(isnan(tempneighs)) = [];
    borderinds = setdiff(tempneighs,find(temp));
    thissizes = zeros(size(thismap));
    thissizes(borderinds) = 1;
end
    
    
    



if ~strcmp(tokens{1},prev)
maps2(:,end+1) = 0;
sizes2(:,end+1) = 0;
col = col+1;
prev = tokens{1};
end
maps2(:,col) = maps2(:,col) + thismap;
sizes2(:,col) = sizes2(:,col) + thissizes;
end
ft_write_cifti_mod([mapname '_' interestingfeature '.dscalar.nii'],maps)
ft_write_cifti_mod([sizesname '_' interestingfeature '.dscalar.nii'],sizes)
maps.data = maps2; maps.mapname = networksfound;
ft_write_cifti_mod([mapname '_' interestingfeature '_combined.dscalar.nii'],maps)
sizes.data = sizes2; sizes.mapname = maps.mapname;
ft_write_cifti_mod([sizesname '_' interestingfeature '_combined.dscalar.nii'],sizes)