function consensusmap = consensus_infomap_results_func(assignmentsfile,ciftitemplatefile,mincol)

mat = load(assignmentsfile);
outputmat = mat;
%outputmat(:,end) = mat(:,end);

final_numcolors = zeros(size(mat,2),1);
final_numunassigned = zeros(size(mat,2),1);

final_numcolors(end) = nnz(unique(outputmat(:,end))>0);
final_numunassigned(end) = nnz(outputmat(:,end)==-1);

for col = [size(mat,2)-1 : -1 : 1]

    col1colors = unique(mat(:,col)); col1colors(col1colors<1) = [];
    col2colors = unique(outputmat(:,col+1)); col2colors(col2colors<1) = [];
    col2sizes = zeros(length(col2colors),1);
    for i=1:length(col2colors)
        col2sizes(i) = nnz(outputmat(:,col+1)==col2colors(i));
    end
    
    [ign sorti] = sort(col2sizes,1,'descend');
    sortcols2 = col2colors(sorti);
    col1colorsavailable = col1colors;
    
    if length(col2colors) > length(col1colors)
        1;
    end
    
    for i = 1:length(sortcols2)
        
        if ~isempty(col1colorsavailable)
        
        col1sizes = zeros(length(col1colorsavailable),1);
        for j = 1:length(col1colorsavailable)
            col1sizes(j) = nnz(mat(outputmat(:,col+1)==sortcols2(i),col)==col1colorsavailable(j));
        end
        [maxval maxi] = max(col1sizes);
        if maxval>0
            outputmat(mat(:,col)==col1colorsavailable(maxi),col) = sortcols2(i);
            col1colorsavailable(maxi) = [];
        else
            if nnz(outputmat(:,col+2:end)==sortcols2(i))==0 
                outputmat(outputmat(:,col+1)==sortcols2(i),col+1) = -1;
            end
        end
        end
    end
    if col==4
        1;
    end
    maxval = max(max(outputmat(:,col+1:end)));%max(col2colors);
    for i = 1:length(col1colorsavailable)
        maxval = maxval+1;
        outputmat(mat(:,col)==col1colorsavailable(i),col) = maxval;
    end
    
final_numcolors(col) = nnz(unique(outputmat(:,col))>0);
final_numunassigned(col) = nnz(outputmat(:,col)==-1);

end

dotindices = strfind(assignmentsfile,'.');
outputfile = [assignmentsfile(1:dotindices(end)-1) '_regularized.txt'];
dlmwrite(outputfile,outputmat,' ');
            

power_surf_colormap = [.67 .67 .67;.67 .67 .67;1 0 0;0 0 .6;1 1 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;1 1 .8;0 .4 0;.25 .25 .25];
if max(unique(outputmat)) < (size(power_surf_colormap,1)-2)
    power_surf_colormap_touse = power_surf_colormap(1:max(unique(outputmat))+2,:);
else
    power_surf_colormap_touse = [power_surf_colormap ; repmat(power_surf_colormap(end,:),(max(unique(outputmat)) - (size(power_surf_colormap,1)-2)),1)];
end

 imagesc(sortrows(outputmat,[size(outputmat,2):-1:1]))
 colormap(power_surf_colormap_touse)
 title('Regularized Assignments')
%figure
%plot([1:size(mat,2)],final_numcolors,'r-')
%figure
%plot([1:size(mat,2)],final_numunassigned,'b-')

% Adjust for non-hierarchical features

allcolors = unique(outputmat(:,mincol:end)); allcolors(allcolors<1)=[];
maxcolor = max(unique(outputmat));
for col = mincol:(size(outputmat,2)-1)
    for color = allcolors'
        if (nnz(outputmat(:,col+1)==color) > 0) && (length(intersect(find(outputmat(:,col)==color),find(outputmat(:,col+1)==color))) < (.6*length(find(outputmat(:,col)==color)))) %&& ((nnz(outputmat(:,col)==color) / size(outputmat,1)) > .01)
            matchcolors = unique(outputmat(outputmat(:,col)==color,col+1));
            matchcolors(matchcolors==color)=[]; matchcolors(matchcolors<1)=[];
            
            mainmatchsize = length(intersect(find(outputmat(:,col)==color),find(outputmat(:,col+1)==color)));
            
            matchcolorsizes = zeros(length(matchcolors),1);
            for matchcolornum = 1:length(matchcolors)
                matchcolor = matchcolors(matchcolornum);
                matchcolorsizes(matchcolornum) = nnz(outputmat(outputmat(:,col)==color,col+1)==matchcolor);
            end
            [secondmatchsize maxi] = max(matchcolorsizes);
            
            if (mainmatchsize + secondmatchsize) > .9*(length(find(outputmat(:,col)==color)))
            
            biggestmatchcolor = matchcolors(maxi);
            maxcolor = maxcolor+1;
            outputmat(logical((outputmat(:,col)==color).*(outputmat(:,col+1)==biggestmatchcolor)),1:col) = maxcolor;
            end
        end
    end
end

% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments

consensusmap = outputmat(:,mincol);


unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = outputmat(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end


% Clean up tiny pieces

minsize = 10;


bufsize = 524288;
ciftiheadertext = textread(ciftitemplatefile,'%s','delimiter','\r','bufsize',bufsize);

for row = 1:length(ciftiheadertext)
    thisline = ciftiheadertext{row};
    
    if ~isempty(strfind(thisline,'<BrainModel IndexOffset="0" IndexCount="'))
        indexloc = strfind(thisline,'IndexCount="');
        quoteloc = strfind(thisline,'"'); quoteloc(quoteloc< (indexloc+length('IndexCount="'))) = []; nextquoteloc = quoteloc(1);
        numverts = str2num(thisline(indexloc+length('IndexCount="') : nextquoteloc-1));
        if numverts < 30000
            iscifti = 1;
        elseif numverts > 30000
            iscifti = 2;
        end
    
    if ~isempty(strfind(thisline,'ModelType="CIFTI_MODEL_TYPE_SURFACE" BrainStructure="CIFTI_STRUCTURE_CORTEX_RIGHT"'))
        hem = 'R';
    elseif ~isempty(strfind(thisline,'ModelType="CIFTI_MODEL_TYPE_SURFACE" BrainStructure="CIFTI_STRUCTURE_CORTEX_LEFT"'))
        hem = 'L';
    end
    
    break
    end
end

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    
end

giftimap = zeros(length(mask),1);
giftimap(logical(mask)) = consensusmap(1:nnz(mask));

allcolors = unique(giftimap(giftimap > 0));

for color = allcolors'
    clusteredmetric = zeros(size(giftimap));
    thiscolorverts = find(giftimap==color);
    for vertex = thiscolorverts'
        
        %find the neighbors of this vertex
        vertexneighbors = neighbors(vertex,:);
        
        %find which of those neighbors also pass the thresholds
        vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
        
        %find if those neighbors have already been assigned different cluster values
        uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
        uniqueneighborvals(uniqueneighborvals==0) = [];
        
        %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
        if isempty(uniqueneighborvals)
            clusteredmetric(vertexneighbors_thiscolor) = vertex;
            %if there is only one previous cluster identifier present, make all the neighbors that value
        elseif length(uniqueneighborvals)==1
            clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
            %if there are multiple cluster identifier values in the neighborhood, merge them into one
        else
            for valuenum = 2:length(uniqueneighborvals)
                clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
            end
        end
    end
    uniqueclustervals = unique(clusteredmetric);
    uniqueclustervals(uniqueclustervals==0) = [];
    
    for clusternum = uniqueclustervals'
        if nnz(clusteredmetric==clusternum) < minsize
            neighborverts = unique(neighbors((clusteredmetric==clusternum),2:7));
            neighborverts(isnan(neighborverts)) = [];
            borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
            borderverts(giftimap(borderverts)<1) = [];
            mode_neighborval = mode(giftimap(borderverts));
            giftimap(clusteredmetric==clusternum) = mode_neighborval;
        end
    end
end

consensusmap(1:nnz(mask)) = giftimap(logical(mask));

