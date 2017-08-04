function edge_distance_comparison(edgemap1name,edgemap2name,threshold,outputfile)

% edgemap1name = 'avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg.func.gii';
% edgemap2name = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg.func.gii';
% threshold = .175;
% outputfile = 'Stevedatavs120_distance.func.gii';

edgemap1 = gifti(edgemap1name);
edgemap1 = edgemap1.cdata;
edgemap1thresh = edgemap1 > threshold;

edgemap2 = gifti(edgemap2name);
edgemap2 = edgemap2.cdata;
edgemap2thresh = edgemap2 > threshold;

edgemap1indices = find(edgemap1thresh);

bufsize=16384;
%caretdir = ['/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K'];
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

outputmetric = ones(size(edgemap2)) .* -1;

for index = edgemap1indices'
    
    distance = 0;
    distancefound = 0;
    nodeneigh = index;
    
    while distancefound == 0
        
        if any(edgemap2thresh(nodeneigh))
            distancefound = 1;
            outputmetric(index) = distance;
        end
        
        newneighs = neighbors(nodeneigh,2:7);
        newneighs = newneighs(1:numel(newneighs));
        nodeneigh = union(nodeneigh, newneighs);
        nodeneigh(isnan(nodeneigh)) = [];
        
        distance = distance+1;
        
    end
end

save(gifti(single(outputmetric)),outputfile);
        
    

