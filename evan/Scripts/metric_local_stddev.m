%function metric_minima(metricname,neighdist,Hem,outputdir,outputname)
tmaskname = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
[subjects tmaskfiles]=textread(tmaskname,'%s%s');

 
 neighdist = 6;
 Hem = 'L';
 %outputdir = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/voxelwiseFC/';
 %outputname = 'smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_AllC_avg_edge_avg_smooth_L_noalone_minima.func.gii';

outputstatsfile = ['/data/cn4/evan/RestingState/FC_Mapping_120/EdgemapQualityMetric.txt'];

delete(outputstatsfile);
fid = fopen(outputstatsfile,'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Subject','StdDev','Minimacount','DataAmount'); %write the output file header
fclose(fid);
dlmwrite(outputstatsfile,' ','-append');
 
medialmaskdata = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/medialmask.func.gii');
medialmaskind = find(~medialmaskdata.cdata);

badnodedata = gifti('/data/cn4/laumannt/fcMapping_redux/all_meanimage_L_32k_fs_LR_750mask.func.gii');
badnodeind = find(badnodedata.cdata);

allbadnodes = union(medialmaskind,badnodeind);

bufsize=16384;
caretdir = '/data/cn4/evan/fsaverage_LR32k/';
%cd(caretdir)
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir 'node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

for s = 1:length(subjects)
    metricname = ['/data/cn4/laumannt/left_hem_edge/' subjects{s} '_avg_edge_avg_smooth_L_noalone_smooth2.55.func.gii'];

metric = gifti(metricname);
allmetric(s,:) = metric.cdata;
tmask{s} = load(tmaskfiles{s});
end

standarddev = zeros(size(allmetric));
minimas = zeros(size(allmetric));

for i = 1:size(allmetric,2)
    
    if ~any(allbadnodes==i);
    
    nodeneigh = i;
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    for n = 1:neighdist
        
        for t = 1:length(curneigh)
            newneigh = [newneigh neighbors(curneigh(t),2:7)];
            newneigh(isnan(newneigh)) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh);
        
    end
    
    [~, badnodeindices, ~] = intersect(nodeneigh,allbadnodes);
    nodeneigh(badnodeindices) = [];
    
    ind = find(nodeneigh==i); %Which nodeneigh is original point
    mins = min(allmetric(:,nodeneigh),[],2);
    %if mini == ind
    
    minimas(:,i) = mins == allmetric(:,i);
    
    standarddev(:,i) = std(allmetric(:,nodeneigh),0,2);
    
    end
    
    
end

standarddev(:,allbadnodes) = [];

for s = 1:length(subjects)

    texttowrite = [subjects{s} '   ' num2str(mean(standarddev(s,:))) '   ' num2str(length(find(minimas(s,:)))) '   ' num2str(length(find(tmask{s})))];
    dlmwrite(outputstatsfile,texttowrite,'-append','delimiter','');
end







