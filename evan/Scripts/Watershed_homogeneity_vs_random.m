hem = 'L';

%colors = [1 0 0; 0 1 0; 0 0 1; 1 1 0; 1 0 1; 0 1 1; .5 .5 .5; .2 .2 .2; .8 .8 .8];
    

gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>750;

threshs = [1.01 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.6 1.7 1.8 1.9 2];

randomPCAs = '/data/cn4/laumannt/parcel_homogeneity/eigval_per_first_rand.mat';
load(randomPCAs);
randomlabels = '/data/cn4/laumannt/parcel_homogeneity/label_rand.mat';
load(randomlabels)
sizes_rand = [];
eigvals_rand = [];
for i = 1:length(eigval_per_first_rand)
    
    for j = 1:100
        
        for k = 1:size(eigval_per_first_rand{i},1)
            
            if nnz((squeeze(single(label_rand(:,j,i)))==k) .* (~gooddata)) == 0
            
                sizes_rand(end+1) = nnz(label_rand(:,j,i)==k);
                eigvals_rand(end+1) = sqrt(FisherTransform(eigval_per_first_rand{i}(k,j)));
            end
        end
    end
end

C_rand = zeros(length(eigvals_rand),3);


%%
allC = [];
allsizes = [];
alleigvals = [];
for threshnum = 1:length(threshs)
    thresh = threshs(threshnum);
    
    colors = [rand rand rand];
    sizes{threshnum} = [];
    eigvals{threshnum} = [];
    C{threshnum} = [];
    
    thiswatershedfilename = ['/data/cn4/evan/Temp/AdjustedWatershedEdges_' num2str(thresh) '_watershedmerge.func.gii'];
    thiswatershed = gifti(thiswatershedfilename); thiswatershed = thiswatershed.cdata;
    thisPCAfilename = ['/data/cn4/evan/Temp/PCA_eigval_per_first_AdjustedWatershedEdges_' num2str(thresh) '_watershedmerge.func.gii'];
    thisPCA = gifti(thisPCAfilename); thisPCA = sqrt(FisherTransform(thisPCA.cdata));
        
    badparcels = unique(thiswatershed .* (~gooddata))';
    for badparcel = badparcels
        thisPCA(thiswatershed==badparcel) = 0;
        thiswatershed(thiswatershed==badparcel) = 0;
    end
    
    for parcel = unique(thiswatershed)'
        if parcel
            allsizes(end+1) = nnz(thiswatershed==parcel);
            alleigvals(end+1) = mean(thisPCA(thiswatershed==parcel));
            
            sizes{threshnum}(end+1) = nnz(thiswatershed==parcel);
            eigvals{threshnum}(end+1) = mean(thisPCA(thiswatershed==parcel));
            C{threshnum}(end+1,:) = colors;
            
            allC(end+1,:) = colors;
        end
    end
    
end
    
    
    