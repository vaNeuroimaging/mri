assignmentsfile = ['rawassn_minsize5.txt'];
parcelfilename = {'/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_L_watershedmerge_0.6_tweaked.func.gii','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_R_watershedmerge_0.6_tweaked.func.gii'};
%{'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_L_nosmooth_watershedmerge_0.4_tweaked.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_R_nosmooth_watershedmerge_0.4_tweaked.func.gii'}; %{'/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_edge_L_watershedmerge_0.45_tweaked.func.gii','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_edge_R_watershedmerge_0.45_tweaked.func.gii'};
templatefile = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';
outputfilestem = 'Poldrome_subsurf_LR_minsize5_consensus';
minsize = 5;

mincol = 10;

recolor = [];
%recolor = [1 10; 2 9; 3 3; 4 1; 5 7; 6 4.5; 7 6; 8 16; 9 5; 10 15; 11 2; 12 2.5; 13 11.5; 14 12;15 17];

%% Regularize

mat = load(assignmentsfile); 
%mat = mat(:,1:32);
outputmat = mat;

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
%             if nnz(outputmat(:,col+2:end)==sortcols2(i))==0 
%                 outputmat(outputmat(:,col+1)==sortcols2(i),col+1) = -1;
%             end
        end
        end
    end

   
    maxval = max(max(outputmat(:,col+1:end)));
    for i = 1:length(col1colorsavailable)
        foundamatch = 0;
        thiscolorinds = find(mat(:,col)==col1colorsavailable(i));
        for highercol = (col+1):size(mat,2)
            thiscol_colors = unique(outputmat(:,highercol)); thiscol_colors(thiscol_colors<1) = [];
            for thiscolor = thiscol_colors(:)'
                numoverlapping = length(intersect(thiscolorinds,find(outputmat(:,highercol)==thiscolor)));
                if ((numoverlapping/length(thiscolorinds)) > .8) && ((numoverlapping/length(find(outputmat(:,highercol)==thiscolor))) > .8)
                    outputmat(thiscolorinds,col) = thiscolor;
                    foundamatch = 1;
                end
            end
            if foundamatch
                break
            end
        end
        if ~foundamatch
            maxval = maxval+1;
            outputmat(thiscolorinds,col) = maxval;
        end
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

figure
 imagesc(sortrows(outputmat,[size(outputmat,2):-1:1]))
 colormap(power_surf_colormap_touse)
 title('Regularized Assignments')
% figure
% plot([1:size(mat,2)],final_numcolors,'r-')
% figure
% plot([1:size(mat,2)],final_numunassigned,'b-')

%% Adjust for non-hierarchical features
% 
% 
% 
% 
% allcolors = unique(outputmat(:,mincol:end)); allcolors(allcolors<1)=[];
% maxcolor = max(unique(outputmat));
% for col = mincol:(size(outputmat,2)-1)
%     for color = allcolors'
%         if (nnz(outputmat(:,col+1)==color) > 0)  && (length(intersect(find(outputmat(:,col)==color),find(outputmat(:,col+1)==color))) < (.3*length(find(outputmat(:,col)==color)))) %&& ((nnz(outputmat(:,col)==color) / size(outputmat,1)) > .01)
%             matchcolors = unique(outputmat(outputmat(:,col)==color,col+1));
%             matchcolors(matchcolors==color)=[]; matchcolors(matchcolors<1)=[];
%             
%             mainmatchsize = length(intersect(find(outputmat(:,col)==color),find(outputmat(:,col+1)==color)));
%             
%             matchcolorsizes = zeros(length(matchcolors),1);
%             for matchcolornum = 1:length(matchcolors)
%                 matchcolor = matchcolors(matchcolornum);
%                 matchcolorsizes(matchcolornum) = nnz(outputmat(outputmat(:,col)==color,col+1)==matchcolor);
%             end
%             [secondmatchsize maxi] = max(matchcolorsizes);
%             
%             if (mainmatchsize + secondmatchsize) > .9*(length(find(outputmat(:,col)==color)))
%             
%             biggestmatchcolor = matchcolors(maxi);
%             maxcolor = maxcolor+1;
%             outputmat(logical((outputmat(:,col)==color).*(outputmat(:,col+1)==biggestmatchcolor)),1:col) = maxcolor;
%             end
%         end
%     end
% end

%% Create initial consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments

consensusmap = outputmat(:,mincol);


unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = outputmat(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end

% for col = mincol:size(outputmat,2)
%     tempconsensus = consensusmap;
%     for parcelnum = 1:length(consensusmap)
%         if (nnz(outputmat(:,col)==outputmat(parcelnum,col)) < (nnz(consensusmap==consensusmap(parcelnum)))/2) && (outputmat(parcelnum,col)~=0)
%             tempconsensus(parcelnum) = outputmat(parcelnum,col);
%         end
%     end
%     consensusmap = tempconsensus;
% end

% 
% %% Check to make sure the selected column isn't badly representitive of the surrounding thresholds
% if mincol>1
%     for node = 1:length(consensusmap)
%         lowervals = outputmat(node,1:(mincol-1)); lowervals(lowervals<1) = [];
%         highervals = outputmat(node,(mincol+1):end); highervals(highervals<1) = [];
%         commonvals = intersect(lowervals,highervals);
%         commonvals(commonvals==consensusmap(node)) = [];
%         thisval_numsurrounding = [];
%         for i = 1:length(commonvals)
%             if (nnz(lowervals==commonvals(i)) > 1 && nnz(highervals==commonvals(i)) > 1)
%                 thisval_numsurrounding(i) = nnz(lowervals==commonvals(i)) + nnz(highervals==commonvals(i));
%             end
%         end
%         [maxnum maxi] = max(thisval_numsurrounding);
%         if ~isempty(maxi)
%             consensusmap(node) = commonvals(maxi);
%         end
%     end
% end
    


%% Reapply size threshold

consensusvals = unique(consensusmap);
for val = consensusvals(:)'
    if nnz(consensusmap==val) < minsize
        toosmallinds = find(consensusmap==val);
        for smallindex = toosmallinds(:)'
            thisassignments = outputmat(unassignedindex,mincol:end);
            thisassignments(thisassignments<1) = [];
            thisassignments(thisassignments==val) = [];
            if ~isempty(thisassignments)
                consensusmap(smallindex) = thisassignments(1);
            end
        end
        
    end
end

%% Move values around to make things prettier

vals_inconsensus = unique(consensusmap);
morevals_in_matrix = setdiff(unique(outputmat),vals_inconsensus);
for val = 1:length(morevals_in_matrix)
    outputmat(outputmat==morevals_in_matrix(val)) = max(vals_inconsensus) + val;
end

for val = 1:max(vals_inconsensus);
    if (nnz(consensusmap==val)==0) && (val<max(vals_inconsensus))
        consensusmap(consensusmap==max(vals_inconsensus)) = val;
        outputmat(outputmat==max(vals_inconsensus)) = val;
        vals_inconsensus = unique(consensusmap);
    end
end


%% Write parcel output



temp_consensusmap = consensusmap;
temp_outputmat = outputmat;
for i = 1:size(recolor,1)
    temp_consensusmap(consensusmap==recolor(i,1)) = recolor(i,2);
    temp_outputmat(outputmat==recolor(i,1)) = recolor(i,2);
end
consensusmap = temp_consensusmap;
outputmat = temp_outputmat;
%outputmat(:,1:mincol-1) = -1;

figure
imagesc(sortrows(outputmat,[size(outputmat,2):-1:1]))
colormap(power_surf_colormap_touse)
title('Recolored regularized assignments')


templatedata = gifti(templatefile);


outputdata = zeros(size(templatedata.cdata,1),1);
crossthresh_outputdata = zeros(size(templatedata.cdata,1),size(outputmat,2));

hems = {'L','R'};
for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
    
    parcels = gifti(parcelfilename{hemnum}); parcels = parcels.cdata(logical(mask));
    
    parcelIDs = unique(parcels(parcels>0));
    if hemnum==1
        nLparcelIDs = length(parcelIDs);
    end
    
    for parcelnum = 1:length(parcelIDs)
        outputdata(find(parcels==parcelIDs(parcelnum)) + (strcmp('R',hem)*ncortexLverts)) = consensusmap(parcelnum + (strcmp('R',hem)*nLparcelIDs));
        for col = 1:size(crossthresh_outputdata,2)
            crossthresh_outputdata(find(parcels==parcelIDs(parcelnum)) + (strcmp('R',hem)*ncortexLverts),col) = outputmat(parcelnum + (strcmp('R',hem)*nLparcelIDs),col);
        end
    end
end
dlmwrite([outputfilestem '.txt'],consensusmap)
cifti_write_wHDR(outputdata,templatefile,outputfilestem)
cifti_write_wHDR(crossthresh_outputdata,templatefile,[outputfilestem '_allthresh'])
