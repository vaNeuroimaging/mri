metricname = '/data/cn4/evan/RestingState/ModifiedVoxelwiseOnsurface/fullmatrix/fullmatrix_Tk0005to0017in0003_S1to1_xd20_BI_INFMAP/ClusterAssignments_kden0.014.func.gii';

metricdata = gifti(metricname);
metricdata = metricdata.cdata;


bufsize=16384;
%caretdir = ['/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K'];
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

border_metric = zeros(size(metricdata));

for i = 1:length(metricdata)
    
    nodeneigh = neighbors(i,2:7);
    nodeneigh(logical(isnan(nodeneigh))) = [];
    neighborvals = metricdata(nodeneigh);
    if length(unique(neighborvals))>1
        border_metric(i) = 1;
    end
end
    
save(gifti(single(border_metric)),[metricname(1:end-9) '_borders.func.gii']);
 
    
            

