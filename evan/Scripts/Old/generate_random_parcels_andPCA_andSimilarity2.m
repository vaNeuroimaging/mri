function generate_random_parcels_andPCA_andSimilarity2(watershedname,iterations,avgcrosscorrname,iscifti,hem,outfileprefix)
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
    
    minimametric(allbadvertices) = 0;
    
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
        waterid = watersheds(waternum);
        if nnz(label==waterid) < minsize
            
            minimametric(label==waterid) = 0;
            
        end
    end
    
    if nnz(minimametric) < numparcels
        error('Not enough parcels left after eliminating small ones!')
    end
    
    while nnz(minimametric) > numparcels
        
        minimaindices = find(minimametric);
        
        randindex = randi(length(minimaindices),1);
        
        minimametric(minimaindices(randindex)) = 0;
        
    end
    
    label = watershed_algorithm_rand(minimametric);
    
    label(allbadvertices) = 0;
    
    multilabel(:,iter) = label;
end

save(gifti(single(multilabel(:,1:10))),[outfileprefix '.func.gii'])

disp(' ')
disp('Evaluating homogeneities of random parcels')
tic
PCAvals = parcel_homogeneity_PCA_multiple_par(multilabel,avgcrosscorr,iscifti,'L');
toc
disp('Evaluating similarities of random parcels')
tic
Simvals = parcel_similarity_multiple_par(multilabel,avgcrosscorr,iscifti,'L');
toc
disp('Calculating averages of random parcels')

for iter = 1:size(multilabel,2)
    
    label = multilabel(:,iter);
    
    watersheds = unique(label);
    watersheds(watersheds==0) = [];
    templabel = zeros(size(label));
    for watershednum = 1:length(watersheds)
        watershed = watersheds(watershednum);
        templabel(label==watershed) = watershednum;
        randomsizes(watershednum,iter) = nnz(label==watershed);
        randomeigvals_per_first(watershednum,iter) = mean(PCAvals(templabel==watershednum,iter));
        randomsimilarities(watershednum,iter) = mean(Simvals(templabel==watershednum,iter));
    end
    
    randommeaneigvals(iter) = mean(randomeigvals_per_first(:,iter));
    randommeansimvals(iter) = mean(randomsimilarities(:,iter));
    
end

disp('Evaluating homogeneities of real parcels')

PCAvals = parcel_homogeneity_PCA_multiple(realwatershed,avgcrosscorr,iscifti,'L');

disp('Evaluating similarities of real parcels')

Simvals = parcel_similarity_multiple(realwatershed,avgcrosscorr,iscifti,'L');

watersheds = unique(realwatershed);
watersheds(watersheds==0) = [];
templabel = zeros(size(realwatershed));
for watershednum = 1:length(watersheds)
    watershed = watersheds(watershednum);
    templabel(realwatershed==watershed) = watershednum;
    realsizes(watershednum) = nnz(realwatershed==watershed);
    realeigvals_per_first(watershednum) = mean(PCAvals(templabel==watershednum));
    realsimilarities(watershednum) = mean(Simvals(templabel==watershednum));
end

realmeaneigval = mean(realeigvals_per_first);
realmeansimval = mean(realsimilarities);

save([outfileprefix '.mat'],'randomsizes', 'randomeigvals_per_first','randomsimilarities','randommeaneigvals','randommeansimvals','realsizes','realeigvals_per_first','realsimilarities','realmeaneigval','realmeansimval')

disp(['Mean difference of real vs random eigvals for ' outfileprefix ': ' num2str(realmeaneigval - mean(randommeaneigvals))])
disp(['Z score of real vs random eigvals for ' outfileprefix ': ' num2str((realmeaneigval - mean(randommeaneigvals))/std(randommeaneigvals))])
disp(['Mean difference of real vs random similarities for ' outfileprefix ': ' num2str(realmeansimval - mean(randommeansimvals))])
disp(['Z score of real vs random similarities for ' outfileprefix ': ' num2str((realmeansimval - mean(randommeansimvals))/std(randommeansimvals))])





%allrandparcels = multilabel;

%end

%save(gifti(single(multilabel)),'/data/cn4/evan/Temp/Temp_randomwater.func.gii');




