hem = 'R';
datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
iscifti = 2;

outnamebase = '/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_vertexwise/DMN_Sal_amPFC_cluster';

 roifilename = 'DMN_Sal_amPFC_ROI.func.gii';
 roi = gifti(roifilename); roi = roi.cdata;

% vertex = 15835;
% roi = zeros(32492,1); roi(vertex) = 1;

networksizeminimum = 10;

datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
iscifti = 2;



if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    roi = roi(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    roi = roi(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    roi = parcels_upsampled;
    ncortexLverts = length(roi);
end


    roiindices = find(roi) + (strcmp(hem,'R') * ncortexLverts);




[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

% evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{1} ' Temp.func.gii']);
% subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
% if ~isempty(tmasklist)
%     tmask = load(tmasks{1});
%     subtimecourse = subtimecourse(:,logical(tmask));
% end
% nsubcortverts = size(subtimecourse,1) - (ncortexLverts+ncortexRverts);


%%

for subnum = 1:length(subjects)

    disp(['Subject ' num2str(subnum)])

    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{subnum});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;

    if subnum==1
        parcelcorrelpatterns = zeros(length(subjects),size(subtimecourse,1));
        nsubcortverts = size(subtimecourse,1) - (ncortexLverts+ncortexRverts);
    end

    parceltimecourses = zeros(1,size(subtimecourse,2));

        parceltimecourses = mean(subtimecourse(roiindices,:),1);
        parcelcorrelpatterns(subnum,:) = paircorr_mod(parceltimecourses',subtimecourse');

end
parcelcorrelpatterns(isnan(parcelcorrelpatterns)) = 0;

output = zeros(32492,1);
output(logical(mask)) = mean(parcelcorrelpatterns(:,(1:nnz(mask)) + (strcmp(hem,'R') * ncortexLverts)),1);
save(gifti(single(output)),[outnamebase '_allsubs_connectivity_' hem '.func.gii'])
%save([outputfolder 'parcelcorrelpatterns_' hem '.mat'],'parcelcorrelpatterns','-v7.3')
%load([outputfolder 'parcelcorrelpatterns_' hem '.mat'])

%%

%parcelnum = find(parcelIDs==parcelID);

%disp(['Parcel ' num2str(parcelnum)])


        Y = pdist(squeeze(parcelcorrelpatterns),'correlation');           % 'pdist' converts the square adjacency matrix to a
        %  1 x n matrix so that the function linkage can construct the tree
        
        clustering = linkage(Y, 'average');      % 'linkage' computes the data to construct the tree
        % 'average' refers to the UPGMA algorithm
        
        clusters = cluster(clustering, 'MaxClust', [1:80]);
        
        subject_correlmat = paircorr_mod(parcelcorrelpatterns');
        
        for numclust = 1:size(clusters,2)
            Qvals(numclust) = M_calc_modularity(clusters(:,numclust),subject_correlmat);
        end
        [maxQval maxQind] = max(Qvals);
        maxQind = 7;
        
        [H,T,perm] = dendrogram(clustering, 0, 'orientation','left', 'colorthreshold', clustering(end-maxQind+1,3)+.000001);  
        orient landscape; 
        
        cophenetic_r = cophenet(clustering, Y);


%%





thisthresh_clusters = clusters(:,maxQind);

communities = unique(thisthresh_clusters);
for comnum = 1:length(communities)
    thisthresh_clusters(thisthresh_clusters==communities(comnum)) = communities(comnum) .* (nnz(thisthresh_clusters==communities(comnum))>=networksizeminimum);
end
communities = unique(thisthresh_clusters);
communities(communities==0) = [];

for comnum = 1:length(communities)
    subindices = (thisthresh_clusters==communities(comnum));
    communitypattern = mean(parcelcorrelpatterns(subindices,:),1);
    
    output = zeros(32492,1);
    output(logical(mask)) = communitypattern((1:nnz(mask)) + (strcmp(hem,'R') * ncortexLverts));
    save(gifti(single(output)),[outnamebase '_cluster' num2str(comnum) '_connectivity_' hem '.func.gii'])
    
end




