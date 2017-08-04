
%metricname = '/data/cn4/evan/RestingState/FC_Mapping_120/smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_AllC_avg_edge_avg_smooth_L_noalone.func.gii';
%'/data/cn4/evan/RestingState/FC_Mapping_120/lOT/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg.func.gii';
neighdist = 3;
mindist = [];
Hem = 'L';
outputdir = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/voxelwiseFC/';

nummaxima = 15;
maximathresh = .3;

outputsinglemetrics = 1;
combinedoutputname = 'minima_allmetrics.func.gii';
metricdir = outputdir;

metrics = dir([metricdir '/Avg_*.func.gii']);



for metricnum = 1:length(metrics)
    
    metricname = metrics(metricnum).name;
    disp(metricname)
    
    outputname = ['minima_' metricname];


bufsize=16384;
caretdir = '/data/cn4/evan/fsaverage_LR32k/';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir 'node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

metric = gifti(metricname);
metric = metric.cdata;

if metricnum == 1
    combinedmaximametric = zeros(size(metric));
end


maximametric = zeros(size(metric));
for i = 1:length(metric)
    
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
    
    ind = find(nodeneigh==i); %Which nodeneigh is original point
    [maxval maxi] = max(metric(nodeneigh));
    %if mini == ind
    if maxval == metric(i)
        maximametric(i) = 1;
    end
    
end

if ~isempty(mindist)
    maximaindices = find(maximametric);
    distances = gifti(['/data/cn4/evan/fsaverage_LR32k/VertexDistances.' Hem '.func.gii']);
    distances = distances.cdata(logical(maximametric),logical(maximametric));
    [xind,yind] = find(distances<mindist);
    
    for tooclosenum = 1:length(xind)
        if metric(maximaindices(xind(tooclosenum))) < metric(maximaindices(yind(tooclosenum)))
            maximametric(maximaindices(xind(tooclosenum))) = 0;
        else
            maximametric(maximaindices(yind(tooclosenum))) = 0;
        end
    end
    
end

if ~isempty(nummaxima)
    maximaindices = find(maximametric);
    maximavals = metric(maximaindices);
    [sortedvals, sortedindices] = sort(maximavals);
    maximametric(maximaindices(sortedindices(1:end-15))) = 0;
end

if ~isempty(maximathresh)
    maximametric(logical(metric<maximathresh)) = 0;
end

if outputsinglemetrics
    save(gifti(single(maximametric)),[outputdir '/' outputname]);
end

combinedmaximametric = combinedmaximametric + maximametric;

end

combinedmaximametric(find(combinedmaximametric)) = 1;

save(gifti(single(combinedmaximametric)),[outputdir '/' combinedoutputname]);





