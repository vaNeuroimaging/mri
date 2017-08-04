tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');


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

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    groupcurv = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.curvature.32k_fs_LR.shape.gii']);
    groupsulc = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sulc.32k_fs_LR.shape.gii']);
    
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])

     mostcommon = gifti(['Variability_' hem '_consensus_dice_mostcommon.func.gii']); mostcommon = mostcommon.cdata;
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
% clusters = zeros(66697,1);
% for ID = IDs(:)'
% IDclusters = metric_cluster_cifti(mostcommon,ID-.5,ID+.5,0);
% for i = 1:size(IDclusters,2)
% clusters = clusters + (IDclusters(:,i) * (max(clusters)+1));
% end
% end
% cifti_write_wHDR(clusters,[],'Variability_LR_consensus_dice_mostcommon_clusters')

groupedges = gifti(['Variability_' hem '_consensus_dice_mostcommon_edges.func.gii']); groupedges = groupedges.cdata;

%cifti_to_gifti('Templatematch_dice_bysubject.dtseries.nii','Templatematch_dice_bysubject')

subjectmatches = gifti(['Templatematch_dice_bysubject_' hem '.func.gii']); subjectmatches = subjectmatches.cdata;
%subjectmatches = gifti(['../Poldrome/Templatematch_dice_Poldrome_' hem '.func.gii']); subjectmatches = subjectmatches.cdata;

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
    
    subcurv = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' hem '.curvature.32k_fs_LR.shape.gii']);
    subsulc = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' hem '.sulc.32k_fs_LR.shape.gii']);
    diffcurv = subcurv.cdata - groupcurv.cdata;
    diffsulc = subsulc.cdata - groupsulc.cdata;

    all_diffcurv(:,s) = diffcurv;
    all_diffsulc(:,s) = diffsulc;
    
end

curvcorrels = zeros(32492,1);
curvps = [];
sulccorrels = zeros(32492,1);
sulcps = [];
for ind = edgeinds(:)'
    nonnans = logical(~isnan(edgemovement(ind,:)));
    [curvcorrels(ind)] = corr(abs(edgemovement(ind,nonnans)'),abs(all_diffcurv(ind,nonnans)'));
    [sulccorrels(ind)] = corr(abs(edgemovement(ind,nonnans)'),abs(all_diffsulc(ind,nonnans)'));
end
save(gifti(single(curvcorrels)),['Edgemovement_vCurvature_correlation_' hem '.func.gii'])
save(gifti(single(sulccorrels)),['Edgemovement_vSulcDepth_correlation_' hem '.func.gii'])


edgemovement = edgemovement(logical(groupedges),:);
all_diffcurv = all_diffcurv(logical(groupedges),:);
all_diffsulc = all_diffsulc(logical(groupedges),:);

if hemnum==1
    edgemovementL = edgemovement;
    all_diffcurv_L = all_diffcurv;
    all_diffsulc_L = all_diffsulc;
    clear all_diffcurv all_diffsulc
else
    edgemovement = [edgemovement;edgemovementL];
    all_diffcurv = [all_diffcurv_L;all_diffcurv];
    all_diffsulc = [all_diffsulc_L;all_diffsulc];
end



end