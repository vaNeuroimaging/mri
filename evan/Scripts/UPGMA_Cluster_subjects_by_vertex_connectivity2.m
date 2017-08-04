hems = {'L','R'};

networktemplatefile = '/data/cn4/evan/Test/Network_templates.dtseries.nii';
system(['wb_command -cifti-convert -to-gifti-ext ' networktemplatefile ' Temp.func.gii']);
networkcorrelpattern = gifti('Temp.func.gii'); networkcorrelpattern=networkcorrelpattern.cdata;
networkIDs = [1:17];

%load /data/cn4/evan/RestingState/Consensus/Network_avgconnectivity_xdist20.mat
%networkIDs = [1 2 3 5 7:16];
%networkIDs = [1:16];
%networkcorrelpattern = networkcorrelpattern(:,networkIDs);

Qvalcutoff = .1;

xdist = 20;

divisions = 80;

variancethresh = .2;

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
    tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
    iscifti = 2;
    %hem = 'R';
    
    %kdenval = .2;
    networksizeminimum = 10;
    
    outputfolder = ['/data/cn4/evan/RestingState/Ind_variability/Clustering/'];
    
    
    bufsize=16384;
    caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    % Read in node neighbor file generated from caret -surface-topology-neighbors
    [neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
        neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
        textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
        'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
    neighbors = neighbors+1;
    
    
    
    
    if iscifti==1
        maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
        maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
    elseif iscifti==2
        maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
        maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
        maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
    end
    
    
    [subjects subdata] = textread(datalist,'%s%s');
    if ~isempty(tmasklist)
        [subjects tmasks] = textread(tmasklist,'%s%s');
    end
    
    load /data/cn4/evan/fsaverage_LR32k/Cifti_eroded_euclidean_distances.mat
    dist = dist((1:nnz(mask)) + (ncortexLverts.*strcmp(hem,'R')),:);
    %dist = (dist>xdist);
    
    %distancefile = ['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'];
    %load(distancefile)
    %geo_distances = geo_distances(logical(mask),logical(mask));
    %geo_distances = (geo_distances > xdist);

    
%     evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{1} ' Temp.func.gii']);
%     subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
%     if ~isempty(tmasklist)
%         tmask = load(tmasks{1});
%         subtimecourse = subtimecourse(:,logical(tmask));
%     end
%     nsubcortverts = size(subtimecourse,1) - (ncortexLverts+ncortexRverts);
    
    cd(outputfolder)
    %%
    
    outputclusterID = zeros(nnz(mask),10);
    outputsize = zeros(nnz(mask),10);
    outputmaxQs = zeros(nnz(mask),1);
    
    divisionlength = ceil(nnz(mask)/divisions);
    for divisionnum = 1:divisions
        if divisionnum == divisions;
            divisionindices = (((divisionnum-1)*divisionlength)+1 : nnz(mask)) + (ncortexLverts.*strcmp(hem,'R'));
        else
            divisionindices = (((divisionnum-1)*divisionlength)+1 : (divisionnum*divisionlength)) + (ncortexLverts.*strcmp(hem,'R'));
        end
        

        
        for subnum = 1:length(subjects)
            
            disp(['Loading subject ' num2str(subnum)])
            
%             if divisionnum==1
%                 evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' Subject' num2str(subnum) '.func.gii']);
%             end
            subtimecourse = gifti(['/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_vertexwise/Subject' num2str(subnum) '.func.gii']); subtimecourse = subtimecourse.cdata;
            tmask = load(tmasks{subnum});
            subtimecourse = subtimecourse(:,logical(tmask))';
            subtimecourse(isnan(subtimecourse)) = 0;
            
            if subnum==1
                subjectcorrelpatterns = zeros(length(subjects),length(divisionindices),size(subtimecourse,2));
            end
            
            subjectcorrelpatterns(subnum,:,:) = paircorr_mod(subtimecourse(:,divisionindices),subtimecourse);
            
        end
        subjectcorrelpatterns(isnan(subjectcorrelpatterns)) = 0;
        for vertexnum = 1:length(divisionindices)
            
            disp(['Vertex ' num2str(divisionindices(vertexnum))])
            
            varianceimage = std(squeeze(subjectcorrelpatterns(:,vertexnum,:)));
            varianceverts = varianceimage > variancethresh;
            
            if nnz(varianceverts) > 200
                Y = pdist(squeeze(subjectcorrelpatterns(:,vertexnum,varianceverts)),'correlation');
                
            clustering = linkage(Y, 'average');
            clusters = cluster(clustering, 'MaxClust', [1:80]);
            
            subject_correlmat = paircorr_mod(squeeze(subjectcorrelpatterns(:,vertexnum,:))');
            
            for numclust = 1:size(clusters,2)
                Qvals(numclust) = M_calc_modularity(clusters(:,numclust),subject_correlmat);
            end
            [maxQval maxQind] = max(Qvals);
            end
            
            
            
            if nnz(varianceverts) > 200 && maxQval >= Qvalcutoff
                
                thisthresh_clusters = clusters(:,maxQind);
                
                communities = unique(thisthresh_clusters);
                for comnum = 1:length(communities)
                    thisthresh_clusters(thisthresh_clusters==communities(comnum)) = communities(comnum) .* (nnz(thisthresh_clusters==communities(comnum))>=networksizeminimum);
                end
                communities = unique(thisthresh_clusters);
                communities(communities==0) = [];
            else
                thisthresh_clusters = ones(size(clusters,1),1);
                communities = unique(thisthresh_clusters);
            end
            
            
            
            communityconnections = zeros(length(communities),1);
            communitysize = zeros(length(communities),1);
            for comnum = 1:length(communities)
                subindices = (thisthresh_clusters==communities(comnum));
                communitypattern = squeeze(mean(subjectcorrelpatterns(subindices,vertexnum,:),1));
                
                for networknum = 1:size(networkcorrelpattern,2)
                    posindices = find(networkcorrelpattern(:,networknum) > 0);
                    %posindices = 1:size(networkcorrelpattern,1);
                    xdist_indices_toavoid = find(dist(vertexnum,:)<xdist);
                    indices = setdiff(posindices,xdist_indices_toavoid)';
                    correlations_with_network_patterns(networknum) = paircorr_mod(communitypattern(indices),networkcorrelpattern(indices,networknum));
                end
                [ign maxindex] = max(correlations_with_network_patterns);
                
                communityconnections(comnum) = networkIDs(maxindex);
                communitysize(comnum) = nnz(thisthresh_clusters==communities(comnum));
                
            end
            
            for comnum = 2:length(communityconnections);
                matches = find(communityconnections(1:(comnum-1))==communityconnections(comnum));
                if ~isempty(matches)
                    communitysize(matches) = communitysize(matches)+communitysize(comnum);
                    communityconnections(comnum) = 0;
                end
            end
            indicestoremove = find(communityconnections==0);
            communityconnections(indicestoremove) = [];
            communitysize(indicestoremove) = [];
            
            [ign sortorder] = sort(communitysize,1,'descend');
            
            outputclusterID(divisionindices(vertexnum) - (ncortexLverts.*strcmp(hem,'R')),1:length(communityconnections)) = communityconnections(sortorder);
            outputsize(divisionindices(vertexnum) - (ncortexLverts.*strcmp(hem,'R')),1:length(communityconnections)) = communitysize(sortorder);
            outputmaxQs(divisionindices(vertexnum) - (ncortexLverts.*strcmp(hem,'R'))) = maxQval;
            
        end
    end
    
    outputdata = zeros(size(mask),size(outputclusterID,2)); outputdata(logical(mask),:) = outputclusterID;
    save(gifti(single(outputdata)),['Subject_clustering_vertexwise_allnetworks_' hem '_minQ_' num2str(Qvalcutoff) '.func.gii'])
    
    outputdata = zeros(size(mask),size(outputsize,2)); outputdata(logical(mask),:) = outputsize;
    save(gifti(single(outputdata)),['Subject_clustering_vertexwise_allnetworks_' hem '_minQ_' num2str(Qvalcutoff) '_clustersize.func.gii'])
    
    outputdata = zeros(size(mask),size(outputmaxQs,2)); outputdata(logical(mask),:) = outputmaxQs;
    save(gifti(single(outputdata)),['Subject_clustering_vertexwise_allnetworks_' hem '_maxQs.func.gii'])
    
    Cluster_and_combine(['Subject_clustering_vertexwise_allnetworks_' hem '_minQ_' num2str(Qvalcutoff) '.func.gii'],hem)
    
end


