function parcel_homogeneity_PCA(parcelfilename,avgcrosscorrname,iscifti,hem)
%parcel_homogeneity_PCA(parcelfilename,avgcrosscorrname,iscifti,hem)

parcels = gifti(parcelfilename); parcels = parcels.cdata;


maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;

    if iscifti == 1
        parcels = parcels(mask==0,:);
        maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
    end

if ~ischar(avgcrosscorrname)
    avgcrosscorr = avgcrosscorrname;
    clear avgcrosscorrname

else
    
    try
        
        avgcrosscorr = gifti(avgcrosscorrname);
        avgcrosscorr = avgcrosscorr.cdata;
        
    catch
        avgcrosscorr = cifti_read(avgcrosscorrname);
    end
    
end



avgcrosscorr = avgcrosscorr([1:nnz(mask==0)] + (strcmp(hem,'R') * ncortexLverts),:);


avgcrosscorr(isnan(avgcrosscorr)) = 0;

index = 1;
for map = 1:size(parcels,2)
    watershed = parcels(:,map);
    waternum = unique(watershed);
    waternum(waternum==0) = [];
    for w = 1:length(waternum)
        brainindices{index} = find(watershed==waternum(w));
        mapnum(index) = map;
        eigval_per_first{index} = 0;
        index = index+1;
    end
end

if iscifti
        outputmetric = zeros(nnz(mask==0),size(parcels,2));
else
        outputmetric = zeros(length(mask),size(parcels,2));
end


%% Calculate PCA, determine explained variance first component
for parcel = 1:length(brainindices)

    
        watercorr = avgcrosscorr(brainindices{parcel},:);
        nanmat = isnan(watercorr);
        nancolsum = sum(nanmat,2);
        watercorr((nancolsum>0),:)= [];
        
        if size(watercorr,1) > 2
            
            [ign ign2 eigvals_per] = PCA_reduction(watercorr','comps',2);
            eigval_per_first{parcel} = eigvals_per(1);
            

        end
    
end


for parcel = 1:length(eigval_per_first)
    outputmetric(brainindices{parcel},mapnum(parcel)) = eigval_per_first{parcel};
end
    
if iscifti
    temp = outputmetric;
    outputmetric = zeros(length(mask),size(outputmetric,2));
    outputmetric(mask==0,:) = temp;
end

slashes = strfind(parcelfilename,'/');
if numel(slashes) > 0
    lastslashloc = slashes(end);
else
    lastslashloc = 0;
end

save(gifti(single(outputmetric)),[parcelfilename(1:lastslashloc) 'PCA_eigval_per_first_' parcelfilename(lastslashloc+1:end)])
end

function [Y V_use eigvals_per] = PCA_reduction(mat,reductype,reducnum)
%TOL

%Calculate covariance of data
covmat = cov(mat);

%Calculate eigenvectors and eigenvalues
[V D] = eig(covmat);

eigvals = diag(D);

[eigvals_sort ind] = sort(eigvals,'descend');

V_sort = V(:,ind);

eigvals_per = (eigvals_sort./sum(eigvals_sort))*100;

switch reductype
    
    case 'comps'
    
        V_use = V_sort(:,1:reducnum);
        
    case 'percentvar'
        
        count = 1;
        persum = 0;
        while persum < reducnum
            
            persum = persum + eigvals_per(count);
            count = count + 1;
        end
        
        V_use = V_sort(:,1:count);
end

%Dimensionality reduction
Y = V_use'*mat';
end
