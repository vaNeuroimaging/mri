function opposed = neighbors_opposed(neighbors)

opposed = zeros(size(neighbors,1),3,2);

for vert = 1:size(neighbors,1)
    level1neighbors = neighbors(vert,2:end);
    level1neighbors(isnan(level1neighbors)) = [];
    
    level1neighbor_pick = level1neighbors(1);
    
    level2neighbors = neighbors(level1neighbor_pick,2:end);
    level2neighbors(isnan(level2neighbors)) = [];
    
    level2neighbors = setdiff(intersect(level2neighbors,level1neighbors),vert);
    
    level2neighbor_1 = level2neighbors(1);
    level3neighbor_forlevel2_1 = neighbors(level2neighbor_1,2:end);
    level3neighbor_forlevel2_1(isnan(level3neighbor_forlevel2_1)) = [];
    level3neighbor_forlevel2_1 = setdiff(intersect(level3neighbor_forlevel2_1,level1neighbors),[vert level1neighbor_pick]);
    
    level2neighbor_2 = level2neighbors(2);
    level3neighbor_forlevel2_2 = neighbors(level2neighbor_2,2:end);
    level3neighbor_forlevel2_2(isnan(level3neighbor_forlevel2_2)) = [];
    level3neighbor_forlevel2_2 = setdiff(intersect(level3neighbor_forlevel2_2,level1neighbors),[vert level1neighbor_pick]);
    
    opposedto_level1neighbor_pick = setdiff(level1neighbors,[level1neighbor_pick level2neighbor_1 level3neighbor_forlevel2_1 level2neighbor_2 level3neighbor_forlevel2_2]);
    
    opposed(vert,1,:) = [level1neighbor_pick opposedto_level1neighbor_pick];
    opposed(vert,2,:) = [level2neighbor_1 level3neighbor_forlevel2_2];
    opposed(vert,3,:) = [level2neighbor_2 level3neighbor_forlevel2_1];
    
end