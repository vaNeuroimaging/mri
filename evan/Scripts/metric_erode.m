function metric_erode(metricname,erosion)

metric = gifti(metricname);
metric = metric.cdata;

bufsize=16384;
%caretdir = ['/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K'];
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

outputmetric = metric;

for column = 1:size(metric,2)
    
    for erosionlevel = 1:erosion
        
        tempmetric = outputmetric(:,column);
    
        indices = find(tempmetric);
        
        for vertex = indices'
            
            nodeneigh = neighbors(vertex,2:7);
            nodeneigh(isnan(nodeneigh)) = [];
            
            if any(tempmetric(nodeneigh) == 0)
                outputmetric(vertex,column) = 0;
            end
        end
    end
end

save(gifti(single(outputmetric)),[metricname(1:end-9) '_eroded' num2str(erosion) '.func.gii']);