function Distance_from_borders(giftispace_infomap,edges,thresh,hem,outputname)

load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])

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

edgesthresh = single(edges > thresh);

netborders = zeros(size(giftispace_infomap));

for i = 1:length(netborders)
    indneighs = neighbors(i,2:end); indneighs(isnan(indneighs)) = [];
    neighvals = giftispace_infomap(indneighs); neighvals(neighvals<1) = [];
    if any(neighvals~=giftispace_infomap(i))
       
        
        newneigh = indneighs;
        curneigh = newneigh;
        
        done = 0;
        numverts = 0;
        if edgesthresh(i) > 0;
            done = 1;
        end
        
    
        while done == 0;
            
            numverts = numverts+1;
            if any(edgesthresh(indneighs)>0)
                done = 1;
            end
            
            for t = 1:length(curneigh)
                newneigh = [newneigh neighbors(curneigh(t),2:7)];
                newneigh(isnan(newneigh)) = [];
            end
            curneigh = setdiff(newneigh,indneighs);
            indneighs = union(indneighs,newneigh);
            
        end
        
        netborders(i) = numverts;
        
    else
        netborders(i) = -1;
    end
end

save(gifti(single(netborders)),outputname);