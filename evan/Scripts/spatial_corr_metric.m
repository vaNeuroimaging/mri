function corr_metric = spatial_corr_metric(metric1,metric2,neighdist,outputdir,outputname)

bufsize=16384;
%caretdir = ['/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K'];
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
cd(caretdir)
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread('node_neighbors.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

corr_metric = zeros(size(metric1));

for i = 1:length(metric1)
    
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
    
    metric1_val = metric1(nodeneigh);
    metric2_val = metric2(nodeneigh);
    corr_metric(i) = paircorr_mod(metric1_val,metric2_val);
    %corr_metric(i) = sum(abs(metric1_val-metric2_val));
    %corr_metric(i) = eta_calc_TL_simple(metric1_val,metric2_val);
end

  %save(gifti(single(corr_metric)),[outputdir '/' outputname '.func.gii']);
 
    
            

