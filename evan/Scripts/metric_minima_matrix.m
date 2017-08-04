function minimamatrix = metric_minima_matrix(matrix,neighdist,hem)

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

%metric(1:12) = 100000;

matrix(mask,:) = 100000;

minimamatrix = zeros(size(matrix));

for i = 1:size(matrix,1)
    
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
    
    nodeneigh(nodeneigh==i) = [];
    nodeneigh = [i nodeneigh];
    
    %ind = find(nodeneigh==i); %Which nodeneigh is original point
    
    [minvals mininds] = min(matrix(nodeneigh,:),[],1);
    
    minimamatrix(i,(mininds==1)) = 1;
    
end



