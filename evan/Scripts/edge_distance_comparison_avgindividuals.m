%function edge_distance_comparison(edgemap1name,edgemap2name,threshold,outputfile)

cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_datalist.txt';
cohortdatafolder = '/data/cn4/laumannt/left_hem_edge/';
hem = 'L';
edgemap1name = '/data/cn4/evan/RestingState/FC_Mapping_120/120_smoothed_edges_L.func.gii';
threshold = .14;
outputfile = '/data/cn4/evan/RestingState/FC_Mapping_120/individual_avg_edgedistance_to_group_edges_L.func.gii';



[ign ign2 ign3 subjects] = textread(cohortfile,'%s %s %s %s');

bufsize=16384;
%caretdir = ['/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K'];
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

edgemap1 = gifti(edgemap1name);
edgemap1 = edgemap1.cdata;
edgemap1thresh = edgemap1 > threshold;
edgemap1indices = find(edgemap1thresh);

outputmetric = ones(size(edgemap1)) .* -1;
outputmetric(edgemap1indices) = 0;

for subject = 1:length(subjects)
    
    string{subject} = ['Running subject ' num2str(subject) ' out of ' num2str(length(subjects))];
    if subject==1; fprintf('%s',string{subject}); else fprintf([repmat('\b',1,length(string{subject-1})) '%s'],string{subject}); end
    
    edgemap2name = [cohortdatafolder '/' subjects{subject} '_avg_edge_avg_smooth_' hem '_noalone_smooth2.55.func.gii'];
    
    edgemap2 = gifti(edgemap2name);
    edgemap2 = edgemap2.cdata;
    edgemap2thresh = edgemap2 > threshold;
    
    distancemetric = zeros(size(edgemap1));
    
    for index = edgemap1indices'
        
        distance = 0;
        distancefound = 0;
        nodeneigh = index;
        
        while distancefound == 0
            
            if any(edgemap2thresh(nodeneigh))
                distancefound = 1;
                distancemetric(index) = distance;
            end
            
            newneighs = neighbors(nodeneigh,2:7);
            newneighs = newneighs(1:numel(newneighs));
            nodeneigh = union(nodeneigh, newneighs);
            nodeneigh(isnan(nodeneigh)) = [];
            
            distance = distance+1;
            
        end
    end
    
    outputmetric = outputmetric + (distancemetric/length(subjects));
    
end

save(gifti(single(outputmetric)),outputfile);



