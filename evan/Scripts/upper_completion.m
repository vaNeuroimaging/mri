function [fUC fUC_smooth plats] = upper_completion(edgemetric)     
% TOL 02/27/13

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
cd(caretdir)
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread('node_neighbors.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

distance = zeros(size(edgemetric));
neighvals = zeros(size(neighbors));

%Avoid NaN
neighvals(1:12,1:6) = edgemetric(neighbors(1:12,1:6));
neighvals(13:end,:) = edgemetric(neighbors(13:end,:));
neighvals(1:12,7) = max(edgemetric)+1;

% Find plateau regions
equalmat = (repmat(edgemetric,[1 size(neighvals,2)])-neighvals)==0;
equalneighinds = sum(equalmat,2)>1;
equalneighs = neighbors(logical(equalneighinds),:);
platnum = 1;
plats = {};
for n = 1:size(equalneighs,1)
    disp(['checking node #' num2str(n) ' of ' num2str(size(equalneighs,1))])
    compareneighs = neighbors(equalneighs(n),logical(equalmat(equalneighs(n),:)));

    platbelong = 0;
    
    for p = 1:length(plats)
        
        %Check if neighbors are members of already made plateau
        if nnz(ismember(compareneighs,plats{p}))>0
            plats{p} = [plats{p} equalneighs(n)];
            platbelong = 1;
            break
        end
    end
    
    if platbelong == 0
        plats{platnum} = equalneighs(n);
        platnum = platnum + 1;
    end
    
end

% Calculate distance to vertex with higher edge value
distance = zeros(size(equalneighs,1),1);
nexthighneighval = zeros(size(equalneighs,1),1);
for n = 1:size(equalneighs,1)
    nodeneighs = equalneighs(n,1);
    
    distval = 1;
    while (distance(n) == 0) && (length(nodeneighs)<(length(edgemetric)-1))
        nodeneighs = neighbors(nodeneighs,2:7);
        nodeneighs(isnan(nodeneighs(:))) = [];
        nodeneighs = unique(nodeneighs(:));
        comp = edgemetric(equalneighs(n,1))<edgemetric(nodeneighs);
        if nnz(comp(:))>0
            distance(n) = distval;
            nexthighneighval(n) = mean(edgemetric(logical(comp)));
        end
        distval = distval + 1;
    end
    disp(['checking node #' num2str(n) ' distance = ' num2str(distance(n))])
end

% Replace edge value on range between plateau value and next highest edge value
% weighted by distance to higher vertex
fUC = edgemetric;
[sortedge sortedgeind] = sort(edgemetric);
for n = 1:size(equalneighs,1)
    if edgemetric(equalneighs(n,1))==max(edgemetric)
    else
        edgediff = (sortedge-edgemetric(equalneighs(n,1)))>0;
        rank = find(edgediff==1);
        nextval = edgemetric(sortedgeind(rank(1)));
        diffval = nextval-edgemetric(equalneighs(n,1));
        %     fUC(equalneighs(n,1)) = edgemetric(equalneighs(n,1))+(1/distance(n))*diffval;
        fUC(equalneighs(n,1)) = edgemetric(equalneighs(n,1))+(1/(distance(n)-(nexthighneighval(n)/max(edgemetric))+1))*diffval;
    end
    
end

% Average node smooth after upper completion
fUC_smooth = fUC;
neighvals(1:12,1:6) = fUC(neighbors(1:12,1:6));
neighvals(13:end,:) = fUC(neighbors(13:end,:));
neighvals(1:12,7) = max(fUC)+1;
fUC_smooth(1:12) = mean(neighvals(1:12,1:6),2);
fUC_smooth(13:end) = mean(neighvals(13:end,:),2);