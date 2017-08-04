%%
hems = {'L','R'};

subthresh = 0;

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    data = gifti(['Variability_' hem '_consensus_dice.func.gii']); data=data.cdata;
    sizes = gifti(['Variability_' hem '_consensus_dice_clustersize.func.gii']);
    threshdata = data .* (sizes.cdata >= subthresh);
    for i = [size(data,2):-1:1]
        if ~any(threshdata(:,i))
            threshdata(:,i) = [];
        end
    end
    save(gifti(single(threshdata)),['Variability_' hem '_consensus_dice_mostcommon.func.gii'])
end

%%
hems = {'L','R'};



for hemnum = 1:length(hems)
    hem = hems{hemnum};

%data = gifti(['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_' hem '.func.gii']);
data = gifti(['Variability_' hem '_consensus_dice_mostcommon.func.gii']);
data = data.cdata(:,1);

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

temp = zeros(size(data));

for i = 1:size(data)

    vertneighs = neighbors(i,2:7); vertneighs(isnan(vertneighs)) = [];
    neighvals = data(vertneighs); %neighvals(neighvals<1) = [];
    if any(neighvals~=data(i)) && data(i) >0%&& (data(i)~=max(neighvals))
        
            temp(i) = data(i);
        
    end
end

save(gifti(single(temp)),['Variability_' hem '_consensus_dice_mostcommon_edges.func.gii'])
%save(gifti(single(temp)),['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_' hem '_edges.func.gii'])
end

%gifti_to_cifti('120_LR_minsize400_recolored_manualconsensus_L_cleaned_edges.func.gii','120_LR_minsize400_recolored_manualconsensus_R_cleaned_edges.func.gii','120_LR_minsize400_recolored_manualconsensus_LR_cleaned_edges')



%%
subthresh = 46;
distance = 8;
sizethreshmm = 100;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

hems = {'L','R'};
for hemnum = 1:length(hems)
    hem = hems{hemnum};
    edges = gifti(['Variability_' hem '_consensus_dice_mostcommon_edges.func.gii']); edges = edges.cdata;
    %edges = gifti(['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_' hem '_edges.func.gii']); edges = edges.cdata;
    variability = gifti(['Variability_' hem '_consensus_dice.func.gii']); variability = variability.cdata(:,2:end);
    numsubs = gifti(['Variability_' hem '_consensus_dice_clustersize.func.gii']); numsubs = numsubs.cdata(:,2:end);
    variability(numsubs<subthresh) = 0;
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])
    
    distvariableregions = zeros(32492,1);
    closevariableregions = ones(32492,1);
    variabilityfromedge = ones(32491,1) .* -1;
    
    IDs = unique(variability); IDs(IDs==0) = [];
    for ID = IDs(:)'
        disp(ID)
        thisID_distvariableregions = zeros(32492,1);
        
        inds = [];
        for col = 1:size(variability,2)
            inds = [inds ;find(variability(:,col)==ID)];
        end
        
        for ind = inds(:)'
            withindist = geo_distances(:,ind) <= distance;
            if ~any(edges(withindist)==ID)
                thisID_distvariableregions(ind) = 1;
            end
        end
        
        
        
        temp = metric_cluster_surfacearea(thisID_distvariableregions,0.5,1.5,sizethreshmm,hem);
        for i = 1:size(temp,2)
            distvariableregions(logical(temp(:,i))) = ID;
        end
        
        
        closevariableregions(logical(thisID_distvariableregions)) = 0;
        
        thisIDedgeinds = find(edges==ID);
        thisIDvariableverts = any(variability==ID,2);
        thisIDvariableverts(logical(thisID_distvariableregions)) = 0;
        thisIDvariableverts_clusters = metric_cluster_surfacearea(thisIDvariableverts,.5,1.5,10,hem);
        clusterborders = zeros(size(thisIDvariableverts_clusters));
        for clusternum = 1:size(thisIDvariableverts_clusters,2)
            
            clusterinds = find(thisIDvariableverts_clusters(:,clusternum)==0);
            for ind = clusterinds(:)'
                indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
                if any(thisIDvariableverts_clusters(indneighs,clusternum)==1) && ~(edges(ind)==ID) && ~(any(edges(indneighs)==ID))
                    clusterborders(ind,clusternum) = 1;
                end
            end
        end
        
        for edgeind = thisIDedgeinds(:)'
            adjacentcluster = 0;
            indneighs = neighbors(edgeind,2:7); indneighs(isnan(indneighs)) = [];
            for clusternum = 1:size(thisIDvariableverts_clusters,2)
                if any(intersect(indneighs,find(thisIDvariableverts_clusters(:,clusternum))))
                    adjacentcluster = clusternum;
                    break
                end
            end
            if adjacentcluster > 0  && (nnz(clusterborders(:,adjacentcluster)) > 0)
                mindist = min(geo_distances(edgeind,logical(clusterborders(:,adjacentcluster))));
                variabilityfromedge(edgeind) = mindist;
            end
        end
            
            
                
        
    end
    
    save(gifti(single(distvariableregions)),['Variabile_regions_' hem '_consensus_dice_distance' num2str(distance) '.func.gii'])
    %save(gifti(single(variabilityfromedge)),['Variabile_regions_' hem '_distance_from_edge.func.gii'])
    
    save(gifti(single(logical(distvariableregions))),['Variable_and_distant_mask_' hem '.func.gii'])
    system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variable_and_distant_mask_' hem '.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Variable_and_distant_mask_' hem '_164.func.gii -largest'])
    distmask164 = gifti(['Variable_and_distant_mask_' hem '_164.func.gii']); distmask164 = distmask164.cdata;
    
    
    save(gifti(single(logical(closevariableregions))),['Variable_and_close_mask_' hem '.func.gii'])
    system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variable_and_close_mask_' hem '.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Variable_and_close_mask_' hem '_164.func.gii -largest'])
    closemask164 = gifti(['Variable_and_close_mask_' hem '_164.func.gii']); closemask164 = closemask164.cdata;
    
    variability = gifti(['Variability_' hem '_consensus_dice_mostcommon.func.gii']);
    save(gifti(single(variability.cdata(:,1))),['Variability_' hem '_consensus_dice_mostcommon.func.gii'])
    system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variability_' hem '_consensus_dice_mostcommon.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Variability_' hem '_consensus_dice_mostcommon_164.func.gii -largest'])
    nonvariable = gifti(['Variability_' hem '_consensus_dice_mostcommon_164.func.gii']); nonvariable = nonvariable.cdata;
    %nonvariable = gifti(['/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_' hem '_164.func.gii']); nonvariable = nonvariable.cdata;
    
    variability164 = gifti(['Variability_' hem '_consensus_diceminsize_' num2str(subthresh) '_combineclusters_164.func.gii']); variability164 = variability164.cdata;
    
    distantvariable = nonvariable; distantvariable(logical(distmask164)) = variability164(logical(distmask164));
    save(gifti(single(distantvariable)),['Variability_' hem '_distantonly_minsize_' num2str(subthresh) '_164.func.gii'])
    
    closevariable = nonvariable; closevariable(logical(closemask164)) = variability164(logical(closemask164));
    save(gifti(single(closevariable)),['Variability_' hem '_closeonly_minsize_' num2str(subthresh) '_164.func.gii'])
    
    
end

gifti_to_cifti(['Variabile_regions_L_consensus_dice_distance' num2str(distance) '.func.gii'],['Variabile_regions_R_consensus_dice_distance' num2str(distance) '.func.gii'],['Variabile_regions_LR_consensus_dice_distance' num2str(distance)])
colored = cifti_read(['Variabile_regions_LR_consensus_dice_distance' num2str(distance) '.dtseries.nii']);
IDs = unique(colored); IDs(IDs==0) = [];
out = zeros(size(colored,1),0);
for ID = IDs(:)'
    temp = metric_cluster_cifti(colored,ID-.5,ID+.5,0);
    out(:,end+1:end+size(temp,2)) = temp;
end

out = out .* repmat(colored,1,size(out,2));
cifti_write_wHDR(out,[],['Variabile_regions_LR_consensus_dice_distance' num2str(distance) 'separated']);
    
    
%%


hems = {'L','R'};
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

segments = zeros(32492*2,0);

gifti_to_cifti('Variability_L_consensus_dice_mostcommon.func.gii','Variability_R_consensus_dice_mostcommon.func.gii','Variability_LR_consensus_dice_mostcommon')
mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon.dtseries.nii');
IDs = unique(mostcommon); IDs(IDs==0) = [];
clusters = zeros(66697,1);
    for ID = IDs(:)'
        IDclusters = metric_cluster_cifti(mostcommon,ID-.5,ID+.5,0);
        for i = 1:size(IDclusters,2)
            clusters = clusters + (IDclusters(:,i) * (max(clusters)+1));
        end
    end
cifti_write_wHDR(clusters,[],'Variability_LR_consensus_dice_mostcommon_clusters')

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])
    
    %mostcommon = gifti(['Variability_' hem '_consensus_dice_mostcommon.func.gii']); mostcommon = mostcommon.cdata;
    %     groupedges = zeros(size(mostcommon));
    %     for i = 1:length(groupedges)
    %         vertneighs = neighbors(i,2:7); vertneighs(isnan(vertneighs)) = [];
    %         neighvals = mostcommon(vertneighs);
    %         if any(neighvals ~= mostcommon(i))
    %             groupedges(i) = mostcommon(i);
    %         end
    %     end
    %     save(gifti(single(groupedges)),['Variability_' hem '_consensus_dice_mostcommon_edges.func.gii']);
    %
    
    
    groupedges = gifti(['Variability_' hem '_consensus_dice_mostcommon_edges.func.gii']); groupedges = groupedges.cdata;
    
    cifti_to_gifti('Templatematch_dice_bysubject.dtseries.nii','Templatematch_dice_bysubject')
    
    subjectmatches = gifti(['Templatematch_dice_bysubject_' hem '.func.gii']); subjectmatches = subjectmatches.cdata;
    
    edgedist = zeros(32492,1);
    edgemovement = zeros(32492,size(subjectmatches,2));
    countofdistsfound = zeros(32492,1);
    
    for s = 1:size(subjectmatches,2)
        disp(num2str(s))
        subdata = subjectmatches(:,s);
        subedges = zeros(size(subdata));
        subedgeneighs = zeros(size(subdata,1),6);
        
        for i = 1:size(subdata,1)
            vertneighs = neighbors(i,2:7); vertneighs(isnan(vertneighs)) = [];
            neighvals = subdata(vertneighs); %neighvals(neighvals<1) = [];
            if any(neighvals~=subdata(i)) && subdata(i) >0%&& (data(i)~=max(neighvals))
                subedges(i) = subdata(i);
                if length(neighvals) == 5
                    neighvals = [neighvals; 0];
                end
                subedgeneighs(i,:) = neighvals;
            end
        end
        
        edgeinds = find(groupedges);
        for ind = edgeinds(:)'
            edgeval = groupedges(ind);
            vertneighs = neighbors(ind,2:7); vertneighs(isnan(vertneighs)) = [];
            neighvals = groupedges(vertneighs); neighvals = unique(neighvals); neighvals(neighvals==0) = []; neighvals(neighvals==edgeval) = [];
            
            if length(neighvals) == 1
                neighvalmatches = any(subedgeneighs==neighvals,2);
                matches = find(neighvalmatches .* (subedges==edgeval));
                [mindist mini] = min(geo_distances(ind,matches));
                matchvert = matches(mini);
                
                if mostcommon(matchvert) == mostcommon(ind)
                    movementdirection = -1;
                elseif mostcommon(matchvert) == neighvals
                    movementdirection = 1;
                else
                    movementdirection = NaN;
                end
                
                
                if mindist < 20
                    countofdistsfound(ind) = countofdistsfound(ind)+1;
                    edgedist(ind) = edgedist(ind) + mindist;
                    
                    edgemovement(ind,s) = mindist * movementdirection;
                    
                else
                    edgemovement(ind,s) = NaN;
                end
                
            end
            
        end
        
    end
    
    edgedist = edgedist ./ countofdistsfound;
    
    save(gifti(single(edgedist)),['Variabile_regions_' hem '_distance_from_edge_bysub.func.gii']);
    
    
    
    
    IDs = unique(groupedges); IDs(IDs==0) = [];
    
    
    
    for edgeID = IDs(:)'
        
        thisID_edges_byneighbor = zeros(32492,max(IDs));
        
        thisID_edgeinds = find(groupedges==edgeID);
        
        for ind = thisID_edgeinds(:)'
            indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
            bordervals = unique(groupedges(indneighs)); bordervals(bordervals==0) = []; bordervals(bordervals==edgeID) = [];
            
            thisID_edges_byneighbor(ind,bordervals) = 1;
        end
        
        movement_meanofsegment{edgeID} = zeros(0,size(subjectmatches,2));
        segmentindex{edgeID} = [];
        
        for neighbornum = 1:size(thisID_edges_byneighbor,2);
            
            thisclustering = metric_cluster(thisID_edges_byneighbor(:,neighbornum),.5,1.5,5);
            
            for numcluster = 1:size(thisclustering,2)
                temp = thisclustering(:,numcluster);
                bigtemp = zeros(32492*2,1);
                bigtemp((32492*(hemnum-1)) + 1 : 32492 + (32492*(hemnum-1))) = temp;
                segments(:,end+1) = bigtemp;
                
                movement_meanofsegment{edgeID}(end+1,:) = nanmean(edgemovement(logical(thisclustering(:,numcluster)),:),1);
                
                segmentindex{edgeID}(end+1) = size(segments,2);
            end
        end
    end
    
    
    %save(gifti(single(segments)),['EdgeSegments_' hem '.func.gii'])
    
    if hemnum==1
        movement_meanofsegment_L = movement_meanofsegment;
        segmentindexL = segmentindex;
    else
        for i = 1:length(movement_meanofsegment)
            movement_meanofsegment{i} = [movement_meanofsegment_L{i} ; movement_meanofsegment{i}];
            segmentindex{i} = [segmentindexL{i} segmentindex{i}];
        end
    end
    
end

gifti_to_cifti('Variabile_regions_L_distance_from_edge_bysub.func.gii','Variabile_regions_R_distance_from_edge_bysub.func.gii','Variabile_regions_LR_distance_from_edge_bysub')

save('Segmentmovementdata.mat','movement_meanofsegment','segmentindex')

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);

segments = segments([(maskL.cdata==0);(maskR.cdata==0)],:);
subcort = zeros(66697 - nnz([(maskL.cdata==0);(maskR.cdata==0)]) , size(segments,2));
segments = [segments;subcort];
cifti_write_wHDR(segments,[],'EdgeSegments_LR');

%%

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');
for i = 1:length(tmasks)
    tmask = load(tmasks{i});
    ntimepoints(i) = nnz(tmask);
end

controlfor_ntimepoints = 0;

make_correlfigs = 0;

clusters = cifti_read('Variability_LR_consensus_dice_mostcommon_clusters.dtseries.nii');
segments = cifti_read('EdgeSegments_LR.dtseries.nii');
load('Segmentmovementdata.mat')

networks = cifti_read('Variability_LR_consensus_dice_mostcommon.dtseries.nii');
IDs = unique(networks); IDs(IDs==0) = [];

bilatinds = [];
unilatinds = [];
numtotaltests = 0;
for edgeID = IDs(:)'
    numtotaltests = numtotaltests + nnz(triu(ones(size(movement_meanofsegment{edgeID},1)),1));
end

poscorrelatededges = zeros(size(segments,1),0);
negcorrelatededges = zeros(size(segments,1),0);
poscorrelatededges_together = zeros(size(segments,2),0);
allps = [];
for runthrough = 1:2
for edgeID = IDs(:)'
    
    
    correls = [];
    ps = [];
    indices = zeros(0,2);
    correlmat = zeros(size(movement_meanofsegment{edgeID},1));
    pmat = zeros(size(movement_meanofsegment{edgeID},1));
    clusterIDs = zeros(size(movement_meanofsegment{edgeID},1),1);
    hemIDs = zeros(size(movement_meanofsegment{edgeID},1),1);
    for i = 1:size(movement_meanofsegment{edgeID},1)
        for j = 1:size(movement_meanofsegment{edgeID},1)
            if i > j
                inds = logical((~isnan(movement_meanofsegment{edgeID}(i,:))) .* (~isnan(movement_meanofsegment{edgeID}(j,:))));
                
                if controlfor_ntimepoints
                    [thesecorrels theseps] = partialcorr(movement_meanofsegment{edgeID}(i,inds)',movement_meanofsegment{edgeID}(j,inds)',ntimepoints(inds)');
                else
                    [thesecorrels theseps] = corrcoef(movement_meanofsegment{edgeID}(i,inds)',movement_meanofsegment{edgeID}(j,inds)');
                    thesecorrels = thesecorrels(1,2);
                    theseps = theseps(1,2);
                end
                correls(end+1) = thesecorrels;
                correlmat(i,j) = thesecorrels;
                ps(end+1) = theseps;
                pmat(i,j) = theseps;
                indices(end+1,:) = [i,j];
                
            end
        end
        clusterIDs(i) = mean(clusters(logical(segments(:,segmentindex{edgeID}(i)))));
        hemIDs(i) = all(find(segments(:,segmentindex{edgeID}(i))) > 29696) + 1;
        
    end
    
    if runthrough==1
        allps = [allps; pmat(logical(tril(ones(size(pmat)),1)))];
    else
        disp(['Network ' num2str(edgeID)])
    %thresh = .05 / numtotaltests;
    %
    %thresh = .05 / nnz(triu(ones(size(pmat)),1));
    %thresh = .05;
    if make_correlfigs
        figure
        
        correlmat = correlmat .* (pmat < thresh);
        [sortedclusterIDs sorti] = sort(clusterIDs,'ascend');
        imagesc(correlmat(sorti,sorti),[-.4 .4]);
        colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
        hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
        combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
        colormap(combined);
        colorbar
        transitions = find(sortedclusterIDs(1:end-1) ~= sortedclusterIDs(2:end));
        hline(transitions+.5,'g')
        vline(transitions+.5,'g')
        
        sortedhemIDs = hemIDs(sorti);
        hemtransitions = find(sortedhemIDs(1:end-1) ~= sortedhemIDs(2:end));
        hline(hemtransitions+.5,'y')
        vline(hemtransitions+.5,'y')
        
        title(['Network ' num2str(edgeID)]);
    end
    
    inds_thatsurvive = find(ps  < thresh);
    
    withincomparisonindices_thatsurvive = indices(inds_thatsurvive,:);
    allindices_thatsurvive = zeros(size(withincomparisonindices_thatsurvive));
    for i = 1:numel(allindices_thatsurvive)
        allindices_thatsurvive(i) = segmentindex{edgeID}(withincomparisonindices_thatsurvive(i));
    end
    rs_thatsurvive = correls(inds_thatsurvive)';
    disp('Relationships that survive correction:')
    disp([allindices_thatsurvive rs_thatsurvive]);
    
%     for i = 1:size(allindices_thatsurvive,1)
%         if (nnz(allindices_thatsurvive(i,:)<=192) == 1) && (nnz(allindices_thatsurvive(i,:)>=193) == 1)
%             bilatinds = [bilatinds allindices_thatsurvive(i,:)];
%         else
%             unilatinds = [unilatinds; allindices_thatsurvive(i,:)];
%         end
%     end
    
    %crossID_allindices_thatsurvive{edgeID} = allindices_thatsurvive;
    
    
    
    for i = 1:size(allindices_thatsurvive,1)
        if rs_thatsurvive(i) > 0
            thiscorrelated = zeros(length(IDs),1); thiscorrelated(allindices_thatsurvive(i,:)) = 1;
            matchexisting = 0;
            for existingnum = 1:size(poscorrelatededges_together,2)
                if any(thiscorrelated==poscorrelatededges_together(:,existingnum))
                    matchexisting = matchexisting+1;
                    if matchexisting==1
                        poscorrelatededges_together(:,existingnum) = single(logical(poscorrelatededges_together(:,existingnum) + thiscorrelated));
                        firstmatch = existingnum;
                    elseif matchexisting==2
                        poscorrelatededges_together(:,firstmatch) = single(logical(poscorrelatededges_together(:,existingnum) + poscorrelatededges_together(:,firstmatch)));
                        poscorrelatededges_together(:,existingnum) = [];
                    end
                end
            end
            if matchexisting==0
                poscorrelatededges_together(:,end+1) = thiscorrelated;
            end
        end
    end
    
    
    
    
    for i = 1:size(allindices_thatsurvive,1)
        if rs_thatsurvive(i) > 0 ;
            poscorrelatededges(:,end+1) = single(logical(sum(segments(:,allindices_thatsurvive(i,:)),2))) .* networks;
        else
            negcorrelatededges(:,end+1) = single(logical(sum(segments(:,allindices_thatsurvive(i,:)),2))) .* networks;
        end
    end
    
    end
    
end
[thresh ign] = FDR( allps, .05);
end

out = zeros(size(poscorrelatededges,1),size(poscorrelatededges_together,2));
for i = 1:size(out,2)
    out(:,i) = single(logical(sum(segments(:,poscorrelatededges_together(:,i)),2))) .* networks;
end
cifti_write_wHDR(out,[],'EdgeSegments_Poscorrelated_FDRcorrected_together')

cifti_write_wHDR(poscorrelatededges,[],'EdgeSegments_Poscorrelated_FDRcorrected')
cifti_write_wHDR(negcorrelatededges,[],'EdgeSegments_Negcorrelated_FDRcorrected')

% %gifti_to_cifti('Variability_L_consensus_dice_mostcommon_edges.func.gii','Variability_R_consensus_dice_mostcommon_edges.func.gii','Variability_LR_consensus_dice_mostcommon_edges')
% edges = cifti_read('Variability_LR_consensus_dice_mostcommon_edges.dtseries.nii');
%
% %load('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat')
% segments = cifti_read('EdgeSegments_LR.dtseries.nii');
%
% bilatsegments = logical(sum(segments(:,bilatinds),2));
% cifti_write_wHDR((edges .* bilatsegments),[],'EdgeSegments_wBilatcorrelations')
%
% unilatsegments = logical(sum(segments(:,unilatinds),2));
% cifti_write_wHDR((edges .* unilatsegments),[],'EdgeSegments_wUnilatcorrelations')


% unilat_adjacentinds = [];
% unilat_nonadjacentinds = [];
%
% for i = 1:size(unilatinds,1)
%     if mean(edges(logical(segments(:,unilatinds(i,1))))) == mean(edges(logical(segments(:,unilatinds(i,2)))))
%         seg1neighs = unique(neighbors(logical(segments(:,unilatinds(i,1))),2:7)); seg1neighs(isnan(seg1neighs)) = [];
%         seg2neighs = unique(neighbors(logical(segments(:,unilatinds(i,2))),2:7)); seg2neighs(isnan(seg2neighs)) = [];
%         if any(intersect(seg1neighs,seg2neighs))
%             unilat_adjacentinds = [unilat_adjacentinds ; unilatinds(i,:)];
%         else
%             unilat_nonadjacentinds = [unilat_nonadjacentinds; unilatinds(i,:)];
%         end
%     end
% end
%
% cifti_write_wHDR((edges .* logical(sum(segments(:,unilat_adjacentinds),2))),[],'EdgeSegments_wUnilat_adjacent_correlations')
% cifti_write_wHDR((edges .* logical(sum(segments(:,unilat_nonadjacentinds),2))),[],'EdgeSegments_wUnilat_nonadjacent_correlations')

            
%% Variable region connectivity t-tests


percentregion_thresh = .25;

all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8_selected.dtseries.nii');
variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon.dtseries.nii');

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

ncortverts = 59412;

%main = zeros(ncortverts,size(variable_regions,2));
%alternate = zeros(ncortverts,size(variable_regions,2));
alternatecount = zeros(1,size(variable_regions,2));

for s = 1:length(subjects)
    disp(subjects{s})
    tmask = load(tmasks{s});
    subdata = cifti_read(ciftifiles{s});
    subdata = subdata(:,logical(tmask));
    
    for r = 1:size(variable_regions,2)
        
        if s==1
            alternate{r} = zeros(ncortverts,0);
            main{r} = zeros(ncortverts,0);
        end
        
        mainID = mode(mostcommon(logical(variable_regions(:,r)))); 
        altID = mean(all_variable_regions(logical(variable_regions(:,r))));
        if (nnz(networkconnection_bysub(logical(variable_regions(:,r)),s)==altID) / nnz(variable_regions(:,r))) >= percentregion_thresh
            alternatecount(r) = alternatecount(r)+1;
            
            indices = logical(variable_regions(:,r) .* (networkconnection_bysub(:,s)==altID));
            
            corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
            corrpattern(isnan(corrpattern)) = 0;
            
            %alternate(:,r) = alternate(:,r) + corrpattern(1:ncortverts);
            alternate{r}(:,end+1) = corrpattern(1:ncortverts);
            
        else
            
            indices = logical(variable_regions(:,r) .* (networkconnection_bysub(:,s)==mainID));
            
            corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
            corrpattern(isnan(corrpattern)) = 0;
            
            %main(:,r) = main(:,r) + corrpattern(1:ncortverts);
            main{r}(:,end+1) = corrpattern(1:ncortverts);
        end
    end
end

mainout = zeros(size(mostcommon,1),size(variable_regions,2));
altout = zeros(size(mostcommon,1),size(variable_regions,2));
tout = zeros(size(mostcommon,1),size(variable_regions,2));
tthresh = tinv(1-(.05 / ncortverts / size(variable_regions,2)),length(subjects));
for r = 1:size(variable_regions,2)
    
    %     alternate(:,r) = alternate(:,r) ./ alternatecount(r);
    %     main(:,r) = main(:,r) ./ (length(subjects) - alternatecount(r));
    
    mainout(1:ncortverts,r) = mean(main{r},2);
    altout(1:ncortverts,r) = mean(alternate{r},2);
    
    [H,P,CI,STATS] = ttest2(main{r}',alternate{r}');
    tout(1:ncortverts,r) = STATS.tstat;
    
end

disp(alternatecount)

cifti_write_wHDR(mainout,[],'MeanConnectivity_mainID_byregion')
cifti_write_wHDR(altout,[],'MeanConnectivity_alternateID_byregion')
cifti_write_wHDR(tout,[],['MeanConnectivity_T_MainVsAlternateID_byregion_correctedT' num2str(tthresh)])
cifti_write_wHDR((mainout - altout),[],'MeanConnectivity_Diff_MainVsAlternateID_byregion')


%% Vaeiable region characteristics

distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
SA = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas.dtseries.nii');
coords = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_ASCIIformat_coords.dtseries.nii');

Variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');

for i = 1:size(Variable_regions,2)
    inds = find(Variable_regions(:,i));
    disp(['Region ' num2str(i)])
    this_SA = sum(SA(inds));
    disp(['Surface area: ' num2str(this_SA)])

    inddistances = sum(distances(inds,inds),1);
    [ign mini] = min(inddistances);
    centroid = inds(mini);
    disp(['Centroid: ' num2str(coords(centroid,:))])
end

    