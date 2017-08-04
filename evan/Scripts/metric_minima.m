function metric_minima(metricname,neighdist,outputdir,outputname,hem,mindist)
%metric_minima(metricname,neighdist,outputdir,outputname,hem,mindist)

% metricname = '/data/cn4/evan/RestingState/FC_Mapping_120/smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_AllC_avg_edge_avg_smooth_L_noalone.func.gii';
% neighdist = 3;
% Hem = 'L';
% outputdir = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/voxelwiseFC/';
% outputname = 'smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_AllC_avg_edge_avg_smooth_L_noalone_minima.func.gii';


maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = logical(mask.cdata);
    

bufsize=16384;
%caretdir = ['/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K'];
caretdir = '/data/cn4/evan/fsaverage_LR32k/';
%cd(caretdir)
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir 'node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

metric = gifti(metricname);
metric = metric.cdata;

metric(1:12) = 100000;

metric(mask) = 100000;

%[temp metric] = upper_completion(metric);
clear temp;

%save(gifti(single(metric)),[outputdir '/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_uc.func.gii']);


minimametric = zeros(size(metric));
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
    [minval mini] = min(metric(nodeneigh));
    minindices = find(metric(nodeneigh)==minval);
    
    
    if minval == metric(i)
        minimametric(i) = 1;
        if length(minindices) > 1
            
            minindices(logical(minindices==ind)) = [];
            metric(nodeneigh(minindices)) = minval+.00001;

        end
    end
    
end

if exist('mindist')
    minimaindices = find(minimametric);
    
    for i = 1:length(minimaindices)
        nodeneigh = minimaindices(i);
        newneigh = nodeneigh;
        curneigh = newneigh;
    
        for n = 1:mindist
            
            for t = 1:length(curneigh)
                newneigh = [newneigh neighbors(curneigh(t),2:7)];
                newneigh(isnan(newneigh)) = [];
            end
            curneigh = setdiff(newneigh,nodeneigh);
            nodeneigh = union(nodeneigh,newneigh);
            
        end
        
        neighstooclose = find(minimametric(nodeneigh));
        
        if ~isempty(neighstooclose);
            
            nodestooclose = union(minimaindices(i),nodeneigh(neighstooclose));
            [ign, minidx] = min(metric(nodestooclose));
            nodestocut = nodestooclose;
            nodestocut(minidx) = [];
            
            minimametric(nodestocut) = 0;
            
        end
    end
end
        


save(gifti(single(minimametric)),[outputdir '/' outputname]);







