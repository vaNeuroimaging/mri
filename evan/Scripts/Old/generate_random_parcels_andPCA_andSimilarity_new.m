function generate_random_parcels_andPCA_andSimilarity_new(watershedname,iterations,avgcrosscorr,iscifti,hem,outfileprefix)
%generate_random_parcels_andPCA_andSimilarity_new(watershedname,iterations,avgcrosscorrname,iscifti,hem,outfileprefix)

%  watershedname = '120_L_crossthresh_watershedmerge.func.gii';
%  iterations = 200;
%  avgcrosscorr = avgcorr.cdata;
%  iscifti = 2;
%  hem = 'L';
%  outfileprefix = '120_L_crossthresh_watershedmerge_randommatched_withsim';

% thresh = 1.7;
% watershedname = ['Poldrome_wateredge_L_' num2str(thresh) '_watershedmerge.func.gii'];
% iterations = 200;
% avgcrosscorr = avgcorr;
% iscifti = 0;
% hem = 'L';
% outfileprefix = ['Poldrome_wateredge_L_' num2str(thresh) '_watershedmerge_randommatched'];

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

avgcrosscorr(isnan(avgcrosscorr)) = 0;

neighdist = 3;
minsize = 15;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;

erodedmask = gifti(['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']); erodedmask = erodedmask.cdata;

cross_datatype_indices = [[1:32492]' zeros(32492,2)];

for i=1:size(cross_datatype_indices,1)
    if ~mask(i)
        cross_datatype_indices(i,2) = max(cross_datatype_indices(:,2))+1;
    end
    if erodedmask(i)
        cross_datatype_indices(i,3) = max(cross_datatype_indices(:,3))+1;
    end
end

% maskindices = cross_datatype_indices(:,1); maskindices(logical(mask)) = 100000000;
% [ign, sortedmaskindices] = sort(maskindices);
% cross_datatype_indices(:,2) = sortedmaskindices;
% 
% maskindices = cross_datatype_indices(:,1); maskindices(~logical(erodedmask)) = 100000000;
% [ign, sortedmaskindices] = sort(maskindices);
%cross_datatype_indices(:,3) = sortedmaskindices;


gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>900;

allgoodindices = find(mask==0 .* gooddata);
allbadvertices = find(logical(mask + (~gooddata)));

realwatershed = gifti(watershedname); realwatershed = realwatershed.cdata;

realwatershed(allbadvertices) = 0;

parcelIDs = unique(realwatershed); parcelIDs(parcelIDs==0) = [];
for parcelID = parcelIDs'
    if nnz(realwatershed==parcelID) < minsize
        realwatershed(realwatershed==parcelID) = 0;
    end
end

parcelIDs = unique(realwatershed); parcelIDs(parcelIDs==0) = [];

numparcels = nnz(unique(realwatershed));

for parcelnum = 1:length(parcelIDs)
    parcelneighbors{parcelnum} = [];
    parcelindices{parcelnum} = find(realwatershed==parcelIDs(parcelnum));
    
end


zeroindices = setdiff(find(realwatershed==0),find(mask));
for index = zeroindices'
    
    nodeneigh = index;
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    for n = 1:neighdist
        
        for t = 1:length(curneigh)
            newneigh = [newneigh neighbors(curneigh(t),2:7)];
            newneigh(isnan(newneigh)) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh);
        
    end
    
    parcelsnearby = unique(realwatershed(nodeneigh)); parcelsnearby(parcelsnearby==0) = [];
    
    for parcelnum = 1:length(parcelsnearby)
        parcelindicesnearby(parcelnum) = find(parcelIDs==parcelsnearby(parcelnum));
    end
    for parcelnum = 1:length(parcelsnearby)
        parcelneighbors{parcelindicesnearby(parcelnum)} = union(parcelneighbors{parcelindicesnearby(parcelnum)},parcelindicesnearby);
    end
    clear parcelindicesnearby
end




surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);

[phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));

phi_new_orig = phi+pi;
theta_new_orig = theta+pi/2;

disp(['Running ' num2str(iterations) ' random rotations of ' num2str(numparcels) ' parcels'])

rotated_eigvals = zeros(numparcels,iterations);
rotated_similarities = zeros(numparcels,iterations);
fprintf('%27s',' ')


for parcelnum = 1:numparcels
    parcelID = parcelIDs(parcelnum);
    
    parcelneighbors{parcelnum}(parcelneighbors{parcelnum}==parcelnum) = [];
    
%    if parcelID == 8343
    
    realparcelindices = find(realwatershed==parcelID);
    
    this_parcel_iter_done = 0;
    
    while this_parcel_iter_done < iterations;
        indextorotate = realparcelindices(randi(length(realparcelindices),1));
        targetindex = allgoodindices(randi(length(allgoodindices),1));
        
        rotationok = 1;
        phirotation = phi(targetindex)-phi(indextorotate);
        thetarotation = theta(targetindex)-theta(indextorotate);
        
        rotparcelindices = zeros(length(realparcelindices),1);
        
        for n = 1:length(realparcelindices)
            thisparcelindex = realparcelindices(n);
            
            phi_new = mod(phi_new_orig(thisparcelindex)+phirotation,2*pi) - pi;
            theta_new = mod(theta_new_orig(thisparcelindex)+thetarotation,pi) - pi/2;
            
            test = repmat([phi_new theta_new],[32492 1]);
            diff_test = sum(abs([phi theta] - test),2);
            [val rotparcelindices(n)] = min(diff_test);
        end
        
        rotatedparcel = zeros(size(mask));
        rotatedparcel(rotparcelindices) = 1;
        rotneighvals = rotatedparcel(neighbors(13:end,2:7));
        rotneighvals_top = rotatedparcel(neighbors(1:12,2:6));
        rotverts_toadd = [((sum(rotneighvals_top,2) > 3) .* (rotatedparcel(1:12)==0)) ; ((sum(rotneighvals,2) > 4) .* (rotatedparcel(13:end)==0))];
        rotatedparcel = rotatedparcel + rotverts_toadd;
        rotparcelindices = find(rotatedparcel);
        truerotparcelindices = cross_datatype_indices(rotparcelindices,iscifti+1);
        
        
        
        if ~isempty(intersect(rotparcelindices,allbadvertices))
            rotationok = 0;
        else
            
            
            
            for parcelneighbornum = 1:length(parcelneighbors{parcelnum})
                thisparcelneighbor = parcelneighbors{parcelnum}(parcelneighbornum);
                
                neighbor_rotparcelindices{parcelneighbornum} = zeros(length(parcelindices{thisparcelneighbor}),1);
                
                for n = 1:length(parcelindices{thisparcelneighbor})
                    thisparcelindex = parcelindices{thisparcelneighbor}(n);
                    
                    phi_new = mod(phi_new_orig(thisparcelindex)+phirotation,2*pi) - pi;
                    theta_new = mod(theta_new_orig(thisparcelindex)+thetarotation,pi) - pi/2;
                    
                    test = repmat([phi_new theta_new],[32492 1]);
                    diff_test = sum(abs([phi theta] - test),2);
                    [val neighbor_rotparcelindices{parcelneighbornum}(n)] = min(diff_test);
                end
                
                rotatedparcel = zeros(size(mask));
                rotatedparcel(neighbor_rotparcelindices{parcelneighbornum}) = 1;
                rotneighvals = rotatedparcel(neighbors(13:end,2:7));
                rotneighvals_top = rotatedparcel(neighbors(1:12,2:6));
                rotverts_toadd = [((sum(rotneighvals_top,2) > 3) .* (rotatedparcel(1:12)==0)) ; ((sum(rotneighvals,2) > 4) .* (rotatedparcel(13:end)==0))];
                rotatedparcel = rotatedparcel + rotverts_toadd;
                neighbor_rotparcelindices{parcelneighbornum} = find(rotatedparcel);
                
                if ~isempty(intersect(neighbor_rotparcelindices{parcelneighbornum},allbadvertices))
                    rotationok = 0;
                    break
                end
                
            end
        end
        
        
        if logical(rotationok) && length(rotparcelindices) > 2
            
            this_parcel_iter_done = this_parcel_iter_done+1;
            
            fprintf(repmat('\b',1,27))
            fprintf('%27s',['Parcel ' num2str(parcelnum) ': rotation ' num2str(this_parcel_iter_done)])
            
            watercorr = avgcrosscorr(truerotparcelindices,:);
            
            
            [ign ign2 eigvals_per] = PCA_reduction(watercorr','comps',2);
            rotated_eigvals(parcelnum,this_parcel_iter_done) = eigvals_per(1);
            
            
            rotcorrelpattern = mean(avgcrosscorr(truerotparcelindices,:),1);
            
            neighbor_rotcorrelpattern = zeros(length(parcelneighbors{parcelnum}),size(rotcorrelpattern,2));
            
            for parcelneighbornum = 1:length(parcelneighbors{parcelnum})
                trueneighbor_rotparcelindices = cross_datatype_indices(neighbor_rotparcelindices{parcelneighbornum},iscifti+1);
                neighbor_rotcorrelpattern(parcelneighbornum,:) = mean(avgcrosscorr(trueneighbor_rotparcelindices,:),1);
            end
            
            correlations = paircorr_mod(rotcorrelpattern',neighbor_rotcorrelpattern');
            
            rotated_similarities(parcelnum,this_parcel_iter_done) = mean(correlations);
            
        end
        clear neighbor_rotparcelindices
    end
%    end
end
    
disp(' ')
disp('Evaluating homogeneities of real parcels')

PCAvals = parcel_homogeneity_PCA_multiple_par(realwatershed,avgcrosscorr,iscifti,'L');

save(gifti(single(PCAvals)),['PCA_eigval_per_first_' watershedname])

[maxSimilarity meanSimilarity]= parcel_similarity(realwatershed,avgcrosscorr,iscifti,'L',3);

save(gifti(single(maxSimilarity)),['Max_similarity_' watershedname])
save(gifti(single(meanSimilarity)),['Mean_similarity_' watershedname])

% disp('Evaluating similarities of real parcels')
% 
% Simvals = parcel_similarity_multiple(realwatershed,avgcrosscorr,iscifti,'L');

watersheds = unique(realwatershed);
watersheds(watersheds==0) = [];
templabel = zeros(size(realwatershed));
for watershednum = 1:length(watersheds)
    watershed = watersheds(watershednum);
    templabel(realwatershed==watershed) = watershednum;
    realsizes(watershednum) = nnz(realwatershed==watershed);
    realeigvals_per_first(watershednum) = mean(PCAvals(templabel==watershednum));
    realminrot_eigvals(watershednum) = realeigvals_per_first(watershednum) - mean(rotated_eigvals(watershednum,:));
    %Zrealminrot_eigvals(watershednum) = (realeigvals_per_first(watershednum) - mean(rotated_eigvals(watershednum,:)))/std(rotated_eigvals(watershednum,:));
    
    realsimilarities(watershednum) = mean(meanSimilarity(templabel==watershednum));
    realminrot_similarities(watershednum) = realsimilarities(watershednum) - mean(rotated_similarities(watershednum,:));
end


save([outfileprefix '.mat'],'realsizes','realeigvals_per_first','rotated_eigvals','realminrot_eigvals','realsimilarities','rotated_similarities','realminrot_similarities')
% 
disp(['Mean difference of real vs rotated eigvals for ' outfileprefix ': ' num2str(mean(realminrot_eigvals))])
%disp(['Mean Z score of real vs rotated eigvals for ' outfileprefix ': ' num2str(mean(Zrealminrot_eigvals))])
% disp(['Z score of real vs random eigvals for ' outfileprefix ': ' num2str((realmeaneigval - mean(randommeaneigvals))/std(randommeaneigvals))])
 disp(['Mean difference of real vs random similarities for ' outfileprefix ': ' num2str(mean(realminrot_similarities))])
% disp(['Z score of real vs random similarities for ' outfileprefix ': ' num2str((realmeansimval - mean(randommeansimvals))/std(randommeansimvals))])







