%assignmentsfile = ['temp_minsize200.txt'];
%outputfilename = 'Poldrome_consensustest_L.func.gii';
 assignmentsfile = ['120_R_rawassn_minsize200.txt'];
 outputfilename = '120_consensustest_R.func.gii';

cifti = 1;
hem = 'R';






if cifti==1;
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
    
else
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii'];
end
mask = gifti(maskname);
mask = mask.cdata;
if cifti==1
    mask = ~mask;
end




mat = load(assignmentsfile);
mat = mat(1:nnz(mask),:);
outputmat = mat;
%outputmat(:,end) = mat(:,end);

final_numcolors = zeros(size(mat,2),1);
final_numunassigned = zeros(size(mat,2),1);

final_numcolors(end) = nnz(unique(outputmat(:,end))>0);
final_numunassigned(end) = nnz(outputmat(:,end)==-1);

for col = [size(mat,2)-1 : -1 : 1]
    if col==4
        1;
    end
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
%figure
%plot([1:size(mat,2)],final_numcolors,'r-')
%figure
%plot([1:size(mat,2)],final_numunassigned,'b-')

%%


        output = zeros(32492,1);


mincol = 2;

allcolors = unique(outputmat(:,mincol:end)); allcolors(allcolors<1)=[];
maxcolor = max(unique(outputmat));
for col = mincol:(size(outputmat,2)-1)
    for color = allcolors'
        if (nnz(outputmat(:,col+1)==color) > 0)  && (length(intersect(find(outputmat(:,col)==color),find(outputmat(:,col+1)==color))) < (.6*length(find(outputmat(:,col)==color)))) %&& ((nnz(outputmat(:,col)==color) / size(outputmat,1)) > .01)
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
figure
imagesc(sortrows(outputmat,[size(outputmat,2):-1:1]))
colormap(power_surf_colormap_touse)
%%
consensusmap = outputmat(:,mincol);

% for row = 1:length(consensusmap)
%     if (consensusmap(row) > 0) && (nnz(outputmat(row,mincol:end)==consensusmap(row)) < 2) && (any(any(outputmat(:,mincol+2:end)==consensusmap(row))))
%         consensusmap(row) = outputmat(row,mincol+2);
%     end
% end


unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = outputmat(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end

values = unique(consensusmap(consensusmap>0));
outputmap = consensusmap;
% for i = 1:length(values);
%     outputmap(consensusmap==values(i)) = i;
% end



output(logical(mask)) = outputmap;
save(gifti(single(output)),outputfilename)
