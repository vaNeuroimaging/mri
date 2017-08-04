maps = ft_read_cifti_mod('Cluster_probability_maps_sorted_10mm_40sub_match108.dscalar.nii');
cifti_temp = ft_read_cifti_mod('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');
surfacearea = cifti_temp.data;
templates = ft_read_cifti_mod('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
templates = templates.data;
neighbors = cifti_neighbors('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');

networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'};

maps.data(59413:end,:) = [];
intemplate = zeros(size(maps.data,1),1);
notintemplate = zeros(size(maps.data,1),1);


for col = 1:size(maps.data,2)
    disp(col)
    %outputclusters = cifti_cluster_surfacearea(sizes.data(:,col),-.5,.5,0,neighbors);
    
    
    thismap = maps.data(:,col);
    
    tokens = tokenize(maps.mapname{col},':');
    networkname = tokens{1};
    tokens = tokenize(maps.mapname{col},'=');
    tokens2 = tokenize(tokens{2},';');
    patchsize = str2num(tokens2{1}(1:end-3));
    
    
    [maxval maxi] = max(thismap);
    thispatch = false(size(thismap));
    thispatch(maxi) = 1;
    sortedvals = sort(unique(thismap),'descend');
    for i = 1:length(sortedvals)
        mask = (thismap > sortedvals(i)) & (~thispatch);
        newverts = 1;
        while ~isempty(newverts)
            maskpos = find(mask);
            thissizemap_neighs = unique(neighbors(logical(thispatch),2:7));
            thissizemap_neighs(isnan(thissizemap_neighs)) = [];
        
            newverts = intersect(thissizemap_neighs,maskpos);
        
            thispatch(newverts) = 1;
        
            mask = mask & (~thispatch);
            
        end
        if sum(surfacearea(logical(thispatch))) >= patchsize;
           break
        end
    end
    
    
    networknumber = find(strcmp(networkname,networklabels));
    
    if any(templates(logical(thispatch))==networknumber);
        prev = intemplate;
        intemplate(logical(thispatch)) = networknumber;
        
        intemplate(logical(prev) & logical(thispatch)) = -1;
    else
        prev = notintemplate;
        notintemplate(logical(thispatch)) = networknumber;
        
        notintemplate(logical(prev) & logical(thispatch)) = -1;
    end
end

cifti_temp.data = intemplate;
ft_write_cifti_mod(['Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_combined_intemplate'],cifti_temp)
cifti_temp.data = notintemplate;
ft_write_cifti_mod(['Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_combined_notintemplate'],cifti_temp)