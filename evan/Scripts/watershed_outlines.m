function outputmetric = watershed_outlines(metricname)


bufsize=16384;
%caredir = '/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K';
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

if ischar(metricname)
    metricdata = gifti(metricname); metricdata = metricdata.cdata;
else
    metricdata = metricname;
end

outputmetric = zeros(length(metricdata),1);

for i = 1:length(metricdata)
    
    theseneighbors = neighbors(i,2:7);
    theseneighbors(isnan(theseneighbors)) = [];
    
    if metricdata(i) == 0 && any(metricdata(theseneighbors) > 0)
        
        outputmetric(i) = 1;
    end
    
    
end

%save(gifti(outputmetric),[metricname(1:end-9) '_outlines.func.gii']);