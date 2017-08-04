function generate_random_parcels_andPCA(watershedname,iterations,avgcrosscorrname,iscifti,hem,outfileprefix)
%generate_random_parcels_andPCA(watershedname,iterations,avgcrosscorrname,iscifti,hem,outfileprefix)


minsize = 10;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;

gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>750;

allgoodindices = find(mask==0 .* gooddata);
allbadvertices = logical(mask + (~gooddata));

realwatershed = gifti(watershedname); realwatershed = realwatershed.cdata;

realwatershed(allbadvertices) = 0;

numparcels = nnz(unique(realwatershed));

if ~ischar(avgcrosscorrname)
    avgcrosscorr = avgcrosscorrname;
    clear avgcrosscorrname
elseif iscifti
    avgcrosscorr = cifti_read(avgcrosscorrname);
else
    
    try
        
        avgcrosscorr = gifti(avgcrosscorrname);
        avgcrosscorr = avgcrosscorr.cdata;
        
    catch
        avgcrosscorr = cifti_read(avgcrosscorrname);
    end
    
end



bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
neighbors(isnan(neighbors)) = 670;

%for parcelcounti = 1:length(parcelcounts)
    %numparcels = parcelcounts(parcelcounti);

    disp(['Running ' num2str(iterations) ' iterations of ' num2str(numparcels) ' random parcels'])

multilabel = zeros(length(mask),iterations);

for iter = 1:iterations
    randomnums = randperm(length(allgoodindices));
    minimaindices = allgoodindices(randomnums <= (numparcels * 1.5 + 100));
    minimametric = zeros(size(mask));
    minimametric(minimaindices) = 1;
    
    minimagood = 0;
    while minimagood == 0
        if any(any(minimametric(neighbors(logical(minimametric),2:7))))
            
            badminima = minimametric .* any(minimametric(neighbors(:,2:7)),2);
            
            minimametric(logical(badminima)) = 0;
            
            nonminimaindices = setdiff(allgoodindices,find(minimametric));
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
    
    label(allbadvertices) = 0;
    
    
%     parcels_incrap = unique((~gooddata) .* (mask==0) .* label);
%     for thisparcel = parcels_incrap'
%         label(label==thisparcel) = 0;
%     end
    
    
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
    
    borderindices = intersect(find(label==0),allgoodindices)';
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
                
                borderlengths_bigenough = 0;
                
                for i = 1:size(watershedborders(waternum,:),2)
                    if length(watershedborders(waternum,i)) > 2
                        borderlengths_bigenough = 1;
                    end
                end
                
                if borderlengths_bigenough
                    
                    done = 0;
                    while done ==0
                        watershedtomergewith = adjacentwatersheds{waternum}(randi(length(adjacentwatersheds{waternum})));
                        if length(watershedborders{waternum,watersheds==watershedtomergewith}) > 1
                            
                            done = 1;
                            label(label==waterid) = mean(label(origlabel==watershedtomergewith));
                            label(watershedborders{waternum,watersheds==watershedtomergewith}) = 0;
                        end
                    end
                    
                    
                else
                    
                    watershedtomergewith = adjacentwatersheds{waternum}(randi(length(adjacentwatersheds{waternum})));
                    label(label==waterid) = mean(label(origlabel==watershedtomergewith));
                    label(watershedborders{waternum,watersheds==watershedtomergewith}) = 0;
                end
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
    
    label(allbadvertices) = 0;
    
   PCAvals = parcel_homogeneity_PCA(label,avgcrosscorr,iscifti,'/data/cn4/evan/Temp/Temp.func.gii','L');
   
   watersheds = unique(label);
   watersheds(watersheds==0) = [];
   templabel = zeros(size(label));
   for watershednum = 1:length(watersheds)
       watershed = watersheds(watershednum);
       templabel(label==watershed) = watershednum;
       randomsizes(watershednum,iter) = nnz(label==watershed);
       randomeigvals_per_first(watershednum,iter) = mean(PCAvals(templabel==watershednum));
   end
   
   randommeaneigvals(iter) = mean(randomeigvals_per_first(:,iter));
    
   multilabel(:,iter) = templabel;

end

disp(' ')

PCAvals = parcel_homogeneity_PCA(realwatershed,avgcrosscorr,iscifti,'/data/cn4/evan/Temp/Temp.func.gii','L');

watersheds = unique(realwatershed);
watersheds(watersheds==0) = [];
templabel = zeros(size(realwatershed));
   for watershednum = 1:length(watersheds)
       watershed = watersheds(watershednum);
       templabel(realwatershed==watershed) = watershednum;
       realsizes(watershednum) = nnz(realwatershed==watershed);
       realeigvals_per_first(watershednum) = mean(PCAvals(templabel==watershednum));
   end
   
realmeaneigval = mean(realeigvals_per_first);

save(gifti(single(multilabel)),[outfileprefix '.func.gii'])
save([outfileprefix '.mat'],'randomsizes', 'randomeigvals_per_first','randommeaneigvals','realsizes','realeigvals_per_first','realmeaneigval')

disp(['Z score of real vs random eigvals for ' outfileprefix ': ' num2str((realmeaneigval - mean(randommeaneigvals))/std(randommeaneigvals))])





%allrandparcels = multilabel;

%end

%save(gifti(single(multilabel)),'/data/cn4/evan/Temp/Temp_randomwater.func.gii');
    
    
    
    
