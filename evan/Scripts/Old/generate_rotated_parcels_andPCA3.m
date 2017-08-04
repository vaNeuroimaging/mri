function generate_rotated_parcels_andPCA3(watershedname,iterations,avgcrosscorr,iscifti,hem,outfileprefix)
%generate_rotated_parcels_andPCA(watershedname,iterations,avgcrosscorrname,iscifti,hem,outfileprefix)

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
gooddata = gooddata>750;

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



surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);

[phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));

phi_new_orig = phi+pi;
theta_new_orig = theta+pi/2;

disp(['Running ' num2str(iterations) ' random rotations of ' num2str(numparcels) ' parcels'])

rotated_eigvals = zeros(numparcels,iterations);

%fprintf('%27s',' ')

%


for parcelnum = 1:numparcels
    parcelID = parcelIDs(parcelnum);
    
    realparcelindices = find(realwatershed==parcelID);
    
    bigphi_theta = repmat([phi(realparcelindices) theta(realparcelindices)],[1 1 length(phi)]);
    bigphi_new_orig = repmat(reshape(phi_new_orig,[1 1 length(phi)]),[length(realparcelindices) 1 1]);
    bigtheta_new_orig = repmat(reshape(theta_new_orig,[1 1 length(phi)]),[length(realparcelindices) 1 1]);
    this_parcel_iter_done = 0;
    
    while this_parcel_iter_done < iterations;
        indextorotate = realparcelindices(randi(length(realparcelindices),1));
        targetindex = allgoodindices(randi(length(allgoodindices),1));
        
        phirotation = phi(targetindex)-phi(indextorotate);
        thetarotation = theta(targetindex)-theta(indextorotate);
        
        %rotparcelindices = zeros(length(realparcelindices),1);
       phi_new = mod(bigphi_new_orig+phirotation,2*pi) - pi;
       theta_new = mod(bigtheta_new_orig+thetarotation,pi) - pi/2;
        
        %phi_theta_new = [reshape(phi_new,[1 1 length(phi)]) reshape(theta_new,[1 1 length(phi)])];
        
        
        diffs = sum(abs(bigphi_theta - [phi_new theta_new]),2);
        [ign rotparcelindices] = min(diffs,[],3);
        
        
%         for n = 1:length(realparcelindices)
%             thisparcelindex = realparcelindices(n);
%             
%             phi_new = mod(phi_new_orig(thisparcelindex)+phirotation,2*pi) - pi;
%             theta_new = mod(theta_new_orig(thisparcelindex)+thetarotation,pi) - pi/2;
%             
%             test = repmat([phi_new theta_new],[32492 1]);
%             diff_test = sum(abs([phi theta] - test),2);
%             [val rotparcelindices(n)] = min(diff_test);
%         end
        
        rotatedparcel = zeros(size(mask));
        rotatedparcel(rotparcelindices) = 1;
        rotneighvals = rotatedparcel(neighbors(13:end,2:7));
        rotneighvals_top = rotatedparcel(neighbors(1:12,2:6));
        rotverts_toadd = [((sum(rotneighvals_top,2) > 3) .* (rotatedparcel(1:12)==0)) ; ((sum(rotneighvals,2) > 4) .* (rotatedparcel(13:end)==0))];
        rotatedparcel = rotatedparcel + rotverts_toadd;
        rotparcelindices = find(rotatedparcel);
        
        if isempty(intersect(rotparcelindices,allbadvertices)) && length(rotparcelindices) > 2
            
            truerotparcelindices = cross_datatype_indices(rotparcelindices,iscifti+1);
            
            this_parcel_iter_done = this_parcel_iter_done+1;
            
            string{this_parcel_iter_done} = ['Parcel ' num2str(parcelnum) ': rotation ' num2str(this_parcel_iter_done)];
            if this_parcel_iter_done==1; fprintf('%s',string{this_parcel_iter_done}); else fprintf([repmat('\b',1,length(string{this_parcel_iter_done-1})) '%s'],string{this_parcel_iter_done}); end
            
            watercorr = avgcrosscorr(truerotparcelindices,:);
            
            
            [ign ign2 eigvals_per] = PCA_reduction(watercorr','comps',2);
            rotated_eigvals(parcelnum,this_parcel_iter_done) = eigvals_per(1);
            
            
        end
        
    end
    fprintf(repmat('\b',1,length(string{this_parcel_iter_done})))
%    end
end
    
disp(' ')
disp('Evaluating homogeneities of real parcels')

PCAvals = parcel_homogeneity_PCA_multiple_par(realwatershed,avgcrosscorr,iscifti,hem);

%save(gifti(single(PCAvals)),['PCA_eigval_per_first_' watershedname])

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
    Zrealminrot_eigvals(watershednum) = (realeigvals_per_first(watershednum) - mean(rotated_eigvals(watershednum,:)))/std(rotated_eigvals(watershednum,:));
    %realsimilarities(watershednum) = mean(Simvals(templabel==watershednum));
end


save([outfileprefix '.mat'],'realsizes','realeigvals_per_first','rotated_eigvals','realminrot_eigvals')
% 
disp(['Mean difference of real vs rotated eigvals for ' outfileprefix ': ' num2str(mean(realminrot_eigvals))])
disp(['Mean Z score of real vs rotated eigvals for ' outfileprefix ': ' num2str(mean(Zrealminrot_eigvals))])

plot(realsizes,realeigvals_per_first,'r.',realsizes,mean(rotated_eigvals,2)','b.')
[H,P,CI,STATS] = ttest(realeigvals_per_first,mean(rotated_eigvals,2)');
disp(['Paired t-test of real vs rotated eigvals for ' outfileprefix ': t = ' num2str(STATS.tstat) ', p = ' num2str(P)])

% disp(['Z score of real vs random eigvals for ' outfileprefix ': ' num2str((realmeaneigval - mean(randommeaneigvals))/std(randommeaneigvals))])
% disp(['Mean difference of real vs random similarities for ' outfileprefix ': ' num2str(realmeansimval - mean(randommeansimvals))])
% disp(['Z score of real vs random similarities for ' outfileprefix ': ' num2str((realmeansimval - mean(randommeansimvals))/std(randommeansimvals))])







