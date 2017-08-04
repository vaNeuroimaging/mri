function generate_rotated_parcellation_andPCA(watershedname,iterations,cov_corr,rotations,iscifti,hem,outfileprefix)
%generate_rotated_parcellation_andPCA(watershedname,iterations,cov_corr,rotations,iscifti,hem,outfileprefix)

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
nonanneighbors = neighbors; nonanneighbors(isnan(neighbors)) = 380;

%avgcrosscorr(isnan(avgcrosscorr)) = 0;

minsize = 15;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;

erodedmask = gifti(['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']); erodedmask = erodedmask.cdata;

% divisions = 10;
% vertsperdivision = floor(size(avgcrosscorr,1)/divisions);
% cov_corr = zeros(size(avgcrosscorr,1));
% for i = 1:divisions
%     
%     for j = 1:divisions
%         %if j >= i
%         disp([i j])
%         vertices{j} = ((i-1)*vertsperdivision+1) : (i*vertsperdivision*(i<divisions) + size(avgcrosscorr,1)*(i==divisions));
%     
%         a = bsxfun(@minus, avgcrosscorr(vertices{i},:)', mean(avgcrosscorr(vertices{i},:)'));
%         b = bsxfun(@minus, avgcrosscorr(vertices{j},:)', mean(avgcrosscorr(vertices{j},:)'));
%         
%         mag_a = sqrt(sum(a.^2, 1));
%         mag_b = sqrt(sum(b.^2, 1));
%         
%         cov_corr(vertices{i},vertices{j}) = (a' * b);
%         %end
%     end
% end




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

realwatershed = gifti(watershedname); realwatershed = realwatershed.cdata;

gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>750;

allgoodindices = find(mask==0 .* gooddata);
allbadvertices = find(logical(mask + (~gooddata)));



%realwatershed(allbadvertices) = 0;

parcelIDs = unique(realwatershed); parcelIDs(parcelIDs==0) = [];
for parcelID = parcelIDs'
    if nnz(realwatershed==parcelID) < (nnz(realwatershed(allbadvertices)==parcelID) + minsize)
        realwatershed(realwatershed==parcelID) = 0;
    end
end


disp('Evaluating homogeneities of real parcels')

PCAvals = parcel_homogeneity_PCA_cov(realwatershed,cov_corr,iscifti,hem);
watershedslashloc = strfind(watershedname,'/');
if isempty(watershedslashloc)
    watershedslashloc = 0;
end

save(gifti(single(PCAvals)),['PCA_eigval_per_first_' watershedname(watershedslashloc(end)+1:end)])


parcelIDs = unique(realwatershed); parcelIDs(parcelIDs==0) = [];

numparcels = nnz(unique(realwatershed));



load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])

%geo_distances_frombad = geo_distances(allbadvertices,:);



surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
spherecoords = sphere.vertices;

[phi theta r] = cart2sph(spherecoords(:,1), spherecoords(:,2),spherecoords(:,3));

phi_new_orig = phi+pi;
theta_new_orig = theta+pi/2;

disp(['Running ' num2str(iterations) ' random rotations of ' num2str(numparcels) ' parcels'])

rotated_eigvals = zeros(numparcels,iterations);

%fprintf('%27s',' ')

%

rotmap = zeros(32492,iterations);
test = zeros(32492,3);

for iternum = 1:iterations
    
    if isempty(rotations)
    
        xrot = rand * 2*pi;
        yrot = rand * 2*pi;
        zrot = rand * 2*pi;

        
    else
        
        xrot = rotations.xrot(iternum);
        yrot = rotations.yrot(iternum);
        zrot = rotations.zrot(iternum);
        
    end

rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];

    

for parcelnum = 1:numparcels
    
    
    parcelID = parcelIDs(parcelnum);
    
    realparcelindices = find(realwatershed==parcelID);
    
    
    %find parcel center
    meanX = mean(spherecoords(realparcelindices,1));
    meanY = mean(spherecoords(realparcelindices,2));
    meanZ = mean(spherecoords(realparcelindices,3));
    coord = [meanX meanY meanZ];
    sphere_coords = [spherecoords(realparcelindices,1) spherecoords(realparcelindices,2) spherecoords(realparcelindices,3)];
    rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
    dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
    [y indval] = min(dist_coord);
    centerind = realparcelindices(indval);
    
    maxdist = max(geo_distances(centerind,realparcelindices));
    
%     valid_targets = find(~(sum((geo_distances_frombad<=maxdist),1)));
%     if isempty(valid_targets)
%         valid_targets = find(~(sum((geo_distances_frombad<=(maxdist/2)),1)));
%     end
    
    %         bigphi_theta = repmat([phi(realparcelindices) theta(realparcelindices)],[1 1 length(phi)]);
    %     bigphi_new_orig = repmat(reshape(phi_new_orig,[1 1 length(phi)]),[length(realparcelindices) 1 1]);
    %     bigtheta_new_orig = repmat(reshape(theta_new_orig,[1 1 length(phi)]),[length(realparcelindices) 1 1]);
    
    
        string{parcelnum} = ['Rotation ' num2str(iternum) ': parcel ' num2str(parcelnum)];
        if parcelnum==1; fprintf('%s',string{parcelnum}); else fprintf([repmat('\b',1,length(string{parcelnum-1})) '%s'],string{parcelnum}); end
        %targetindex = valid_targets(randi(length(valid_targets),1));
        
%         phirotation = phi(targetindex)-phi(centerind);
%         thetarotation = theta(targetindex)-theta(centerind);
        
        %rotparcelindices = zeros(length(realparcelindices),1);
        %        phi_new = mod(bigphi_new_orig+phirotation,2*pi) - pi;
        %        theta_new = mod(bigtheta_new_orig+thetarotation,pi) - pi/2;
        %
        %         %phi_theta_new = [reshape(phi_new,[1 1 length(phi)]) reshape(theta_new,[1 1 length(phi)])];
        %
        %
        %         diffs = sum(abs(bigphi_theta - [phi_new theta_new]),2);
        %         [ign rotparcelindices] = min(diffs,[],3);
        
        indexcoords = spherecoords(realparcelindices,:)';
            xrotcoords = rotmat_x * indexcoords;
            xyrotcoords = rotmat_y * xrotcoords;
            xyzrotcoords = rotmat_z * xyrotcoords;
        
        clear rotparcelindices
        for n = 1:length(realparcelindices)
            
            test(:,1) = xyzrotcoords(1,n); test(:,2) = xyzrotcoords(2,n); test(:,3) = xyzrotcoords(3,n);
            
            %test = repmat(xyzrotcoords',[32492 1]);
            diff_test = sum(abs(spherecoords - test),2);
            [val rotparcelindices(n)] = min(diff_test);
            if n==indval
                rotcenterind = rotparcelindices(n);
            end
        end
        
        doneadding = 0;
        while doneadding ==0
        rotatedparcel = zeros(size(mask));
        rotatedparcel(rotparcelindices) = 1;
        rotneighvals = rotatedparcel(neighbors(13:end,2:7));
        rotneighvals_top = rotatedparcel(neighbors(1:12,2:6));
        rotverts_toadd = [((sum(rotneighvals_top,2) > 3) .* (rotatedparcel(1:12)==0)) ; ((sum(rotneighvals,2) > 4) .* (rotatedparcel(13:end)==0))];
        if nnz(rotverts_toadd) == 0
            doneadding = 1;
        else
        rotatedparcel = rotatedparcel + rotverts_toadd;
        rotparcelindices = find(rotatedparcel);
        end
        end
        
        
            
        
        if (numel(intersect(rotparcelindices,find(mask))) == 0) %&& (numel(intersect(rotparcelindices,find(~gooddata))) < 15)
        
        samesize = 0;
        count = 0;
        while samesize==0

            deltaverts = length(rotparcelindices) - length(realparcelindices);
            rotatedparcel = zeros(size(rotatedparcel)); rotatedparcel(rotparcelindices) = 1; rotatedparcel(logical(mask)) = 2;
            if sign(deltaverts) == 1
                
                borderverts = find((rotatedparcel==1) .* any(rotatedparcel(nonanneighbors(:,2:end))==0,2));
                
                if length(borderverts) >= deltaverts
                    rotparcelindices = setdiff(rotparcelindices,borderverts(1:deltaverts));
                    samesize = 1;
                else
                    rotparcelindices = setdiff(rotparcelindices,borderverts);
                end
            elseif sign(deltaverts) == -1
                
                borderverts = find((rotatedparcel==0) .* any(rotatedparcel(nonanneighbors(:,2:end))==1,2));
                
                if length(borderverts) >= abs(deltaverts)
                    rotparcelindices = union(rotparcelindices,borderverts(1:abs(deltaverts)));
                    samesize = 1;
                else
                    rotparcelindices = union(rotparcelindices,borderverts);
                end;
            else
                samesize = 1;
            end
        end
        %rotparcelindices = setdiff(rotparcelindices,allbadvertices);
        
        rotatedparcel = zeros(size(rotatedparcel)); rotatedparcel(rotparcelindices) = 1;
        borderverts = find((rotatedparcel==0) .* any(rotatedparcel(nonanneighbors(:,2:end))==1,2));
        
        if (numel(intersect(rotparcelindices,find(mask))) == 0)  && ((numel(borderverts) / numel(rotparcelindices)) < 3) %&& (numel(intersect(rotparcelindices,find(~gooddata))) < 15)
            
            rotmap(rotparcelindices,iternum) = randi(max(parcelIDs),1);
            
            truerotparcelindices = cross_datatype_indices(rotparcelindices,iscifti+1);
            
            watercovcorr = cov_corr(truerotparcelindices(logical(truerotparcelindices)),truerotparcelindices(logical(truerotparcelindices)));
            
            
            eigvals_per = PCA_reduction_cov_onecomp(watercovcorr);
            rotated_eigvals(parcelnum,iternum) = eigvals_per(1);

        else
            
            rotated_eigvals(parcelnum,iternum) = NaN;
        end  
        
        else
            
            rotated_eigvals(parcelnum,iternum) = NaN;
        end
        
        
    end
    fprintf(repmat('\b',1,length(string{parcelnum})))
    %    end
end
    
disp(' ')


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
    %realminrot_eigvals(watershednum) = realeigvals_per_first(watershednum) - mean(rotated_eigvals(watershednum,:));
    %Zrealminrot_eigvals(watershednum) = (realeigvals_per_first(watershednum) - mean(rotated_eigvals(watershednum,:)))/std(rotated_eigvals(watershednum,:));
    %realsimilarities(watershednum) = mean(Simvals(templabel==watershednum));
end


save([outfileprefix '.mat'],'realsizes','realeigvals_per_first','rotated_eigvals')%,'realminrot_eigvals')

save(gifti(single(rotmap(:,1:20))),[outfileprefix '_rotatedmaps.func.gii'])
%save(gifti(single(rotmap)),[outfileprefix '_rotatedmaps.func.gii'])
% 
%disp(['Mean difference of real vs rotated eigvals for ' outfileprefix ': ' num2str(nanmean(realminrot_eigvals))])
%disp(['Mean Z score of real vs rotated eigvals for ' outfileprefix ': ' num2str(nanmean(Zrealminrot_eigvals))])

plot(realsizes,realeigvals_per_first,'r.',realsizes,nanmean(rotated_eigvals,2)','b.')
[H,P,CI,STATS] = ttest(realeigvals_per_first,nanmean(rotated_eigvals,2)');
disp(['Paired t-test of real vs rotated eigvals for ' outfileprefix ': t = ' num2str(STATS.tstat) ', p = ' num2str(P)])

% disp(['Z score of real vs random eigvals for ' outfileprefix ': ' num2str((realmeaneigval - mean(randommeaneigvals))/std(randommeaneigvals))])
% disp(['Mean difference of real vs random similarities for ' outfileprefix ': ' num2str(realmeansimval - mean(randommeansimvals))])
% disp(['Z score of real vs random similarities for ' outfileprefix ': ' num2str((realmeansimval - mean(randommeansimvals))/std(randommeansimvals))])







