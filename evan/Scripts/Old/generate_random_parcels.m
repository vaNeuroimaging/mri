function allrandparcels = generate_random_parcels(parcelcounts,iterations,hem)

%numparcels = 50;
%iterations = 100;
%hem = 'L';

minsize = 10;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
corticalindices = find(mask==0);

gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>750;


bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
neighbors(isnan(neighbors)) = 670;

for parcelcounti = 1:length(parcelcounts)
    numparcels = parcelcounts(parcelcounti);

    disp(['Running ' num2str(iterations) ' iterations of ' num2str(numparcels) ' random parcels'])

multilabel = zeros(length(mask),iterations);

for iter = 1:iterations
    randomnums = randperm(length(corticalindices));
    minimaindices = corticalindices(randomnums <= (numparcels * 1.5 + 100));
    minimametric = zeros(size(mask));
    minimametric(minimaindices) = 1;
    
    minimagood = 0;
    while minimagood == 0
        if any(any(minimametric(neighbors(logical(minimametric),2:7))))
            
            badminima = minimametric .* any(minimametric(neighbors(:,2:7)),2);
            
            minimametric(logical(badminima)) = 0;
            
            nonminimaindices = setdiff(corticalindices,find(minimametric));
            randomnums = randperm(length(nonminimaindices));
            newminimaindices = nonminimaindices(randomnums <= (nnz(badminima)));
            
            minimametric(newminimaindices) = 1;
        else
            minimagood = 1;
        end
    end
    
    string{iter} = ['   Iteration ' num2str(iter)];
    if iter==1; fprintf('%s',string{iter}); else fprintf([repmat('\b',1,length(string{iter-1})) '%s'],string{iter}); end
    
    label = watershed_algorithm_rand(minimametric);
    
    
    parcels_incrap = unique((~gooddata) .* (mask==0) .* label);
    for thisparcel = parcels_incrap'
        label(label==thisparcel) = 0;
    end
    
    
    origlabel = label;
    
    watersheds = unique(label);
    watersheds(watersheds==0) = [];
    
    %mergesneeded = length(watersheds) - numparcels;
    %mergesdone = 0;
    
    
    for waternum = 1:length(watersheds)
        adjacentwatersheds{waternum} = [];
        for waternum2 = 1:length(watersheds)
            watershedborders{waternum,waternum2} = [];
        end
        
    end
    
    borderindices = intersect(intersect(find(label==0),corticalindices),find(gooddata))';
    for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        
        watershedneighbors = unique(borderneighvals);
        
        for waterneighbornum = 1:length(watershedneighbors)
            
            thiswatershedindex(waterneighbornum) = find(watersheds==watershedneighbors(waterneighbornum));
            
            otherwatershedneighbors = watershedneighbors;
            otherwatershedneighbors(otherwatershedneighbors==watershedneighbors(waterneighbornum)) = [];
            
            adjacentwatersheds{thiswatershedindex(waterneighbornum)} = unique([adjacentwatersheds{thiswatershedindex(waterneighbornum)} otherwatershedneighbors']);
        end
        
        if length(watershedneighbors) == 2 && length(borderneighvals)>2
            
            watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
            watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        end
    end
    
    for waternum = 1:length(watersheds)
        waterid = watersheds(waternum);
        if nnz(label==waterid) < minsize
            
            if isempty(adjacentwatersheds{waternum})
                
                label(label==waterid) = 0;
            else
            
                watershedtomergewith = adjacentwatersheds{waternum}(randi(length(adjacentwatersheds{waternum})));
                label(label==waterid) = mean(label(origlabel==watershedtomergewith));
                label(watershedborders{waternum,watersheds==watershedtomergewith}) = 0;
            end
            %mergesdone = mergesdone+1;
                
        end
    end
    
    if nnz(unique(label)) < numparcels
        error('Not enough parcels left after eliminating small ones!')
    end
    
    while nnz(unique(label)) > numparcels
        labelspresent = unique(label); labelspresent(labelspresent==0) = [];
        labeltomerge = labelspresent(randi(length(labelspresent)));
        
        if ~isempty(adjacentwatersheds{watersheds==labeltomerge})
        
            watershedtomergewith = adjacentwatersheds{watersheds==labeltomerge}(randi(length(adjacentwatersheds{watersheds==labeltomerge})));
            
            label(label==labeltomerge) = mean(label(origlabel==watershedtomergewith));
            label(watershedborders{waternum,watersheds==watershedtomergewith}) = 0;
            %mergesdone = mergesdone+1;
            
        end
    end
    
    for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1;
            label(bordervertex) = unique(borderneighvals);
        end
    end
    
   multilabel(:,iter) = label;

end

disp(' ')



allrandparcels{parcelcounti} = multilabel;

end

%save(gifti(single(multilabel)),'/data/cn4/evan/Temp/Temp_randomwater.func.gii');
    
    
    
    
