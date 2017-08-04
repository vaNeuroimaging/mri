function generate_rotated_parcellation_andPCA_nobaddata_generic(watershedname,iterations,cov_corr,rotations,datainds,hem,outfileprefix,baddata)
%generate_rotated_parcellation_andPCA_nobaddata_generic(watershedname,iterations,cov_corr,rotations,datainds,hem,outfileprefix,[baddataname])



bufsize=16384;
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread('/home/data/scripts/Resources/node_neighbors.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
nonanneighbors = neighbors; nonanneighbors(isnan(neighbors)) = 380;

%Set the number of vertices that must be outside baddata regions
minsize = 15;

mask = ~datainds;



cross_datatype_indices = [[1:length(mask)]' zeros(length(mask),1)];

for i=1:size(cross_datatype_indices,1)
    if ~mask(i)
        cross_datatype_indices(i,2) = max(cross_datatype_indices(:,2))+1;
    end
end

realwatershed = gifti(watershedname); realwatershed = realwatershed.cdata;


gooddata = (mask==0) .* (~baddata);


allgoodindices = find((mask==0) .* (~baddata));
allbadvertices = find(logical(mask + baddata));


parcelIDs = unique(realwatershed); parcelIDs(parcelIDs==0) = [];
for parcelID = parcelIDs'
    if nnz(realwatershed==parcelID) < (nnz(realwatershed(allbadvertices)==parcelID) + minsize)
        realwatershed(realwatershed==parcelID) = 0;
    end
end


disp(['Evaluating homogeneities of parcels in ' watershedname])

PCAvals = parcel_homogeneity_PCA_cov_generic(realwatershed,cov_corr,datainds);
watershedslashloc = strfind(watershedname,'/');
if isempty(watershedslashloc)
    watershedslashloc = 0;
end

save(gifti(single(PCAvals)),['PCA_eigval_per_first_' watershedname(watershedslashloc(end)+1:end)])


% parcelIDs = unique(realwatershed); parcelIDs(parcelIDs==0) = [];
% 
% numparcels = nnz(unique(realwatershed));
% 
% 
% 
% geo_distances = ft_read_cifti_mod(['/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Surface_distances_' hem '.dconn.nii']);
% geo_distances = geo_distances.data;
% 
% 
% 
% surfdir = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/';
% sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.coord.gii']);
% spherecoords = sphere.vertices;
% 
% [phi theta r] = cart2sph(spherecoords(:,1), spherecoords(:,2),spherecoords(:,3));
% 
% phi_new_orig = phi+pi;
% theta_new_orig = theta+pi/2;
% 
% disp(['Running ' num2str(iterations) ' random rotations of ' num2str(numparcels) ' parcels'])
% 
% rotated_eigvals = zeros(numparcels,iterations);
% 
% 
% 
% rotmap = zeros(32492,iterations);
% test = zeros(32492,3);
% 
% for iternum = 1:iterations
%     
%     if isempty(rotations)
%     
%         xrot = rand * 2*pi;
%         yrot = rand * 2*pi;
%         zrot = rand * 2*pi;
% 
%         
%     else
%         
%         xrot = rotations.xrot(iternum);
%         yrot = rotations.yrot(iternum);
%         zrot = rotations.zrot(iternum);
%         
%     end
% 
% rotmat_x = [1 0 0;0 cos(xrot) -sin(xrot); 0 sin(xrot) cos(xrot)];
% rotmat_y = [cos(yrot) 0 sin(yrot); 0 1 0; -sin(yrot) 0 cos(yrot)];
% rotmat_z = [cos(zrot) -sin(zrot) 0; sin(zrot) cos(zrot) 0; 0 0 1];
% 
%     
% 
% for parcelnum = 1:numparcels
%     
%     
%     parcelID = parcelIDs(parcelnum);
%     
%     realparcelindices = find(realwatershed==parcelID);
%     
%     
%     %find parcel center
%     meanX = mean(spherecoords(realparcelindices,1));
%     meanY = mean(spherecoords(realparcelindices,2));
%     meanZ = mean(spherecoords(realparcelindices,3));
%     coord = [meanX meanY meanZ];
%     sphere_coords = [spherecoords(realparcelindices,1) spherecoords(realparcelindices,2) spherecoords(realparcelindices,3)];
%     rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
%     dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
%     [y indval] = min(dist_coord);
%     centerind = realparcelindices(indval);
%     
%     maxdist = max(geo_distances(centerind,realparcelindices));
%     
%     
%     
%         string{parcelnum} = ['Rotation ' num2str(iternum) ': parcel ' num2str(parcelnum)];
%         if parcelnum==1; fprintf('%s',string{parcelnum}); else fprintf([repmat('\b',1,length(string{parcelnum-1})) '%s'],string{parcelnum}); end
% 
%         
%         indexcoords = spherecoords(realparcelindices,:)';
%             xrotcoords = rotmat_x * indexcoords;
%             xyrotcoords = rotmat_y * xrotcoords;
%             xyzrotcoords = rotmat_z * xyrotcoords;
%         
%         clear rotparcelindices
%         for n = 1:length(realparcelindices)
%             
%             test(:,1) = xyzrotcoords(1,n); test(:,2) = xyzrotcoords(2,n); test(:,3) = xyzrotcoords(3,n);
%             diff_test = sum(abs(spherecoords - test),2);
%             [val rotparcelindices(n)] = min(diff_test);
%             if n==indval
%                 rotcenterind = rotparcelindices(n);
%             end
%         end
%         
%         doneadding = 0;
%         while doneadding ==0
%         rotatedparcel = zeros(size(mask));
%         rotatedparcel(rotparcelindices) = 1;
%         rotneighvals = rotatedparcel(neighbors(13:end,2:7));
%         rotneighvals_top = rotatedparcel(neighbors(1:12,2:6));
%         rotverts_toadd = [((sum(rotneighvals_top,2) > 3) .* (rotatedparcel(1:12)==0)) ; ((sum(rotneighvals,2) > 4) .* (rotatedparcel(13:end)==0))];
%         if nnz(rotverts_toadd) == 0
%             doneadding = 1;
%         else
%         rotatedparcel = rotatedparcel + rotverts_toadd;
%         rotparcelindices = find(rotatedparcel);
%         end
%         end
%         
% 
%             
%         
%         if (numel(intersect(rotparcelindices,find(mask))) == 0) && (numel(intersect(rotparcelindices,find(~gooddata))) < 15)
%         
%         samesize = 0;
%         count = 0;
%         while samesize==0
% 
%             deltaverts = length(rotparcelindices) - length(realparcelindices);
%             rotatedparcel = zeros(size(rotatedparcel)); rotatedparcel(rotparcelindices) = 1; rotatedparcel(logical(mask)) = 2;
%             if sign(deltaverts) == 1
%                 
%                 borderverts = find((rotatedparcel==1) .* any(rotatedparcel(nonanneighbors(:,2:end))==0,2));
%                 
%                 if length(borderverts) >= deltaverts
%                     rotparcelindices = setdiff(rotparcelindices,borderverts(1:deltaverts));
%                     samesize = 1;
%                 else
%                     rotparcelindices = setdiff(rotparcelindices,borderverts);
%                 end
%             elseif sign(deltaverts) == -1
%                 
%                 borderverts = find((rotatedparcel==0) .* any(rotatedparcel(nonanneighbors(:,2:end))==1,2));
%                 
%                 if length(borderverts) >= abs(deltaverts)
%                     rotparcelindices = union(rotparcelindices,borderverts(1:abs(deltaverts)));
%                     samesize = 1;
%                 else
%                     rotparcelindices = union(rotparcelindices,borderverts);
%                 end;
%             else
%                 samesize = 1;
%             end
%         end
%         
%         rotatedparcel = zeros(size(rotatedparcel)); rotatedparcel(rotparcelindices) = 1;
%         borderverts = find((rotatedparcel==0) .* any(rotatedparcel(nonanneighbors(:,2:end))==1,2));
%         
%         if (numel(intersect(rotparcelindices,find(mask))) == 0)  && ((numel(borderverts) / numel(rotparcelindices)) < 3) && (numel(intersect(rotparcelindices,find(~gooddata))) < 15)
%             
%             rotmap(rotparcelindices,iternum) = randi(max(parcelIDs),1);
%             
%             truerotparcelindices = cross_datatype_indices(rotparcelindices,2);
%             
%             watercovcorr = cov_corr(truerotparcelindices(logical(truerotparcelindices)),truerotparcelindices(logical(truerotparcelindices)));
%             
%             
%             eigvals_per = PCA_reduction_cov_onecomp(watercovcorr);
%             rotated_eigvals(parcelnum,iternum) = eigvals_per(1);
% 
%         else
%             
%             rotated_eigvals(parcelnum,iternum) = NaN;
%         end  
%         
%         else
%             
%             rotated_eigvals(parcelnum,iternum) = NaN;
%         end
%         
%         
%     end
%     %fprintf(repmat('\b',1,length(string{parcelnum})))
% 
% end
%     
% disp(' ')
% 
% 
% watersheds = unique(realwatershed);
% watersheds(watersheds==0) = [];
% templabel = zeros(size(realwatershed));
% for watershednum = 1:length(watersheds)
%     watershed = watersheds(watershednum);
%     templabel(realwatershed==watershed) = watershednum;
%     realsizes(watershednum) = nnz(realwatershed==watershed);
%     realeigvals_per_first(watershednum) = mean(PCAvals(templabel==watershednum));
%     
% end
% 
% 
% save([outfileprefix '.mat'],'realsizes','realeigvals_per_first','rotated_eigvals');







