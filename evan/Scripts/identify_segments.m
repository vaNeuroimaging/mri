function identify_segments(edgesfilename)

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


edges = gifti(edgesfilename); edges = edges.cdata;

edgeverts = find(edges);
divisionverts = [];
for vert = edgeverts(:)'
    vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs)) = [];
    if nnz(edges(vertneighs)) > 2
        divisionverts(end+1) = vert;
    end
end

segnum = 1;
allsegs = zeros(length(edges),1);
thisvert = divisionverts(1);
allsegs(thisvert,segnum) = segnum;

vertsused = [];

while numel(vertsused) < (numel(edgeverts) - numel(divisionverts))
    
    string = ['verts to classify: ' num2str((numel(edgeverts) - numel(divisionverts) - numel(vertsused)))];
    if ~exist('laststring'); 
        fprintf('%s',string); 
    else
        fprintf([repmat('\b',1,length(laststring)) '%s'],string); 
    end
    laststring = string;
    
    vertneighs = neighbors(thisvert,2:7); vertneighs(isnan(vertneighs)) = [];
    edgevertneighs = vertneighs(logical(edges(vertneighs)));
    newedgevertneighs = setdiff(edgevertneighs,union(vertsused,divisionverts));
    
    if ~isempty(newedgevertneighs)    
        thisvert = newedgevertneighs(1);
        allsegs(thisvert,segnum) = segnum;
    else
        for divisionvert = divisionverts(:)'
            vertneighs = neighbors(divisionvert,2:7); vertneighs(isnan(vertneighs)) = [];
            edgevertneighs = vertneighs(logical(edges(vertneighs)));
            if ~isempty(setdiff(edgevertneighs,union(vertsused,divisionverts)));
                thisvert = divisionvert;
                break
            end
        end
    end
        
    
    
    
    if any(divisionverts==thisvert)
        segnum = segnum+1;
        allsegs(:,segnum) = 0;
        allsegs(thisvert,segnum) = segnum;
    else
        vertsused(end+1) = thisvert;
    end
    
end
disp(' ')
outname = [edgesfilename(1:end-9) '_segments.func.gii'];
save(gifti(single(allsegs)),outname)
        
    
    
    
    

