function label = watershed_algorithm_rand(minimametric)
home = pwd;
bufsize=16384;
%caredir = '/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K';
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

labelnum = nnz(minimametric);
labelpos = find(minimametric==1);
label = zeros(size(minimametric));
for j = 1:labelnum;
label(labelpos(j)) = j;
basin{j} = find(label==j);
end

origlabel = zeros(size(label));
%for i = 1:iter
iternum = 1;
while nnz(origlabel ~= label)>0
   % disp(['iteration #' num2str(iternum)])
    origlabel = label;
    for b = 1:length(basin)
        
        nodeneigh = neighbors(basin{b},:);
        nodeneigh = unique(nodeneigh(:));
        nodeneigh(isnan(nodeneigh)) = [];
        newnodeneigh = setdiff(nodeneigh,basin{b});
        
        newnodeneigh_neigh = neighbors(newnodeneigh,:);
        newnodeneigh_neigh(isnan(newnodeneigh_neigh)) = 670; %This is an index in the medial wall
        newnodeneigh_neighvals = label(newnodeneigh_neigh);
        
        basin{b} = [basin{b}; newnodeneigh(sum((newnodeneigh_neighvals~=b) & (newnodeneigh_neighvals~=0),2)==0)];
        
        
        %basin{b} = [basin{b}; newnodeneigh(sum(newnodeneigh_neighvals>0,2)==1)];
        label(basin{b}) = b;
        
    end
%     for j = 1:labelnum
%         label(basin{j}) = j;
%     end
    iternum = iternum + 1;
end



