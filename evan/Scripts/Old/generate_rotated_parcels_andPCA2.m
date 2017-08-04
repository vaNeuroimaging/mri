%function generate_rotated_parcels_andPCA2(watershedname,iterations,avgcrosscorr,iscifti,hem,outfileprefix)
%generate_rotated_parcels_andPCA(watershedname,iterations,avgcrosscorrname,iscifti,hem,outfileprefix)

thresh = 2.6;
watershedname = ['Poldrome_wateredge_L_' num2str(thresh) '_watershedmerge.func.gii'];
iterations = 200;
avgcrosscorr = avgcorr;
iscifti = 0;
hem = 'L';
['Poldrome_wateredge_L_' num2str(thresh) '_watershedmerge_randommatched']

avgcrosscorr(isnan(avgcrosscorr)) = 0;

minsize = 10;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;

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


disp(['Running ' num2str(iterations) '+ rotations of ' num2str(numparcels) ' random parcels'])


rotated_eigvals = zeros(numparcels,iterations);

parcel_iter_done = zeros(numparcels,1);

rotation = 0;
fprintf('%35s',' ')

    
    phi_new = phi+pi;
    theta_new = theta+pi/2;
    
    %Randomize
    randphi = 2*pi.*repmat(rand([1 1 iterations*5]),32492,1);
    randtheta = pi.*repmat(rand([1 1 iterations*5]),32492,1);
    
    phi_new = mod(repmat(phi_new, [1 1 iterations*5])+randphi,2*pi)-pi;
    theta_new = mod(repmat(theta_new, [1 1 iterations*5])+randtheta,2*pi)-pi/2;
    
    bigphi = repmat(phi, [1 1 iterations*5]);
    bigtheta = repmat(theta, [1 1 iterations*5]);
    
    for n = 1:length(phi_new)
        fprintf(repmat('\b',1,35))
        fprintf('%35s',['Remapping rotations of vertex ' num2str(n)])
        test = repmat([phi_new(n,1,:) theta_new(n,1,:)],[32492 1 1]);
        diff_test = sum(abs([bigphi bigtheta] - test),2);
        [val manyind(n,1,:)] = min(diff_test,[],1);
        
    end
    disp(' ')
fprintf('%27s',' ')
while any(parcel_iter_done<iterations)
    
    %disp(['Rotation ' num2str(rotation)])
    rotation = rotation+1;

    ind = squeeze(manyind(:,1,rotation));
    
    rotatedwatershed = realwatershed(ind);
    
    for parcelnum = 1:numparcels
        parcelID = parcelIDs(parcelnum);
        fprintf(repmat('\b',1,27))
        fprintf('%27s',['Rotation ' num2str(rotation) ', parcel ' num2str(parcelnum)])
        
        if (parcel_iter_done(parcelnum) < iterations) && (isempty(intersect(find(rotatedwatershed==parcelID),allbadvertices)))
            
            parcel_iter_done(parcelnum) = parcel_iter_done(parcelnum)+1;

            watercorr = avgcrosscorr((rotatedwatershed==parcelID),:);
            
            %if size(watercorr,1) > 2
                
                [ign ign2 eigvals_per] = PCA_reduction(watercorr','comps',2);
                rotated_eigvals(parcelnum,parcel_iter_done(parcelnum)) = eigvals_per(1);
                
            %end
            
        end
        
    end
end
    

disp('Evaluating homogeneities of real parcels')

PCAvals = parcel_homogeneity_PCA_multiple(realwatershed,avgcrosscorr,iscifti,'L');

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
% disp(['Z score of real vs random eigvals for ' outfileprefix ': ' num2str((realmeaneigval - mean(randommeaneigvals))/std(randommeaneigvals))])
% disp(['Mean difference of real vs random similarities for ' outfileprefix ': ' num2str(realmeansimval - mean(randommeansimvals))])
% disp(['Z score of real vs random similarities for ' outfileprefix ': ' num2str((realmeansimval - mean(randommeansimvals))/std(randommeansimvals))])







