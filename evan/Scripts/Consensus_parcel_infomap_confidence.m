assignmentsfile = ['rawassn_minsize3.txt'];
confidencefile = 'confidence_minsize3.txt';
second_idfile = 'second_id.txt';
secondid_confidencefile = 'secondid_confidence.txt';
alpha = .75;
addlthresh = .1;


%parcelfilename = ['/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_82_L_watershedmerge_0.45_tweaked.func.gii'];
parcelfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii';
hem = 'L';
outputfilename = '120_L_minsize3_consensus_mincol1.func.gii';
outputfilestem = [outputfilename(1:end-9)];
minsize = 3;

mincol = 1;

recolor = [];
%recolor = [1 4; 2 9;3 11;4 1;5 2;6 5;7 7;8 15;9 16;10 4.5;11 17;12 6;13 8;14 6.5; 15 3;16 10.5; 17 10; 18 6];
recolor = [1 1; 2 10;3 3;4 9;5 2;6 15;7 6.5;8 5;9 6;10 18;11 4;12 7;13 11;14 11.5];

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

%% Create consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments

consensusmap = outputmat(:,mincol);


unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = outputmat(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end

for col = mincol:size(outputmat,2)
    tempconsensus = consensusmap;
    for parcelnum = 1:length(consensusmap)
        if (nnz(outputmat(:,col)==outputmat(parcelnum,col)) < (nnz(consensusmap==consensusmap(parcelnum)))/2) && (outputmat(parcelnum,col)~=0)
            tempconsensus(parcelnum) = outputmat(parcelnum,col);
        end
    end
    consensusmap = tempconsensus;
end


%% Check to make sure the selected column isn't badly representitive of the surrounding thresholds
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

for val = 1:max(vals_inconsensus);
    if (nnz(consensusmap==val)==0) && (val<max(vals_inconsensus))
        consensusmap(consensusmap==max(vals_inconsensus)) = val;
        vals_inconsensus = unique(consensusmap);
    end
end




%% Consensusize second assignments
alt_consensusmap = zeros(length(consensusmap),ceil(1/addlthresh));

second_id = load(second_idfile);
second_id_remapped = zeros(size(second_id));

for col = 1:size(second_id,2)
    values = unique(second_id(:,col)); values(isnan(values)) = []; values(values<1) = [];
    for value = values(:)'
        remappedvalue = mode(outputmat(mat(:,col)==value,col));
        second_id_remapped(second_id(:,col)==value,col) = remappedvalue;
    end
end

%% Build matrix of primary and alternate node communities
        
confidence = load(confidencefile);
mean_confidence = mean(confidence,2);

secondid_confidence = load(secondid_confidencefile);

for node = 1:size(consensusmap,1)
    alt_consensusmap(node,1) = consensusmap(node);
    if (consensusmap(node) > 0) && (mean_confidence(node) < alpha)
        altvaluenum = 1;
        thisnode_second_id_values = unique(second_id_remapped(node,:)); thisnode_second_id_values(isnan(thisnode_second_id_values)) = []; thisnode_second_id_values(thisnode_second_id_values<1) = [];
        thisnode_second_id_values(thisnode_second_id_values==consensusmap(node)) = [];
        for value = thisnode_second_id_values(:)'
            if (mean(secondid_confidence(node,(second_id_remapped(node,:)==value))) > addlthresh) && (nnz(second_id_remapped(node,:)==value) > 1)
                altvaluenum = altvaluenum+1;
                alt_consensusmap(node,altvaluenum) = value;
            end
        end
        if ~any(alt_consensusmap(node,2:end))
            alt_consensusmap(node,2) = -1;
        end
    end
end
            



%% Recolor parcels

temp_consensusmap = consensusmap;
temp_outputmat = outputmat;
temp_alt_consensusmap = alt_consensusmap;
for i = 1:size(recolor,1)
    temp_consensusmap(consensusmap==recolor(i,1)) = recolor(i,2);
    temp_outputmat(outputmat==recolor(i,1)) = recolor(i,2);
    temp_alt_consensusmap(alt_consensusmap==recolor(i,1)) = recolor(i,2);
end
consensusmap = temp_consensusmap;
outputmat = temp_outputmat;
alt_consensusmap = temp_alt_consensusmap;

%% Write parcel output

%outputmat(:,1:mincol-1) = -1;

% figure
% imagesc(sortrows(outputmat,[size(outputmat,2):-1:1]))
% colormap(power_surf_colormap_touse)
% title('Recolored regularized assignments')

parcels = gifti(parcelfilename); parcels = parcels.cdata;
outputgifti = zeros(size(parcels));
parcelIDs = unique(parcels(parcels>0));
for parcelnum = 1:length(parcelIDs)
    outputgifti(parcels==parcelIDs(parcelnum)) = consensusmap(parcelnum);
end
dlmwrite([outputfilestem '.txt'],consensusmap)
save(gifti(single(outputgifti)),outputfilename);


system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample ' parcelfilename ' /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC ' parcelfilename(1:end-9) '_164.func.gii -largest'])
parcels_upsampled = gifti([parcelfilename(1:end-9) '_164.func.gii']); parcels_upsampled = parcels_upsampled.cdata;



outputgifti = zeros(size(parcels_upsampled));
parcelIDs = unique(parcels_upsampled(parcels_upsampled>0));
for parcelnum = 1:length(parcelIDs)
    outputgifti(parcels_upsampled==parcelIDs(parcelnum)) = consensusmap(parcelnum);
end
dlmwrite([outputfilestem '.txt'],consensusmap)
save(gifti(single(outputgifti)),[outputfilestem '_164.func.gii']);




sphere = gifti(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii']);
[phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
thetavals = -pi/2 : 1/70 : pi/2;
surf_withlines = zeros(size(phi));
for i = 1:length(thetavals)-1
    theseindices = intersect(find(theta>thetavals(i)), find(theta>thetavals(i+1)));
    surf_withlines(theseindices) = i;
end
surf_withlines = surf_withlines+1;
    

outputgifti = zeros(size(parcels_upsampled));
parcelIDs = unique(parcels_upsampled(parcels_upsampled>0));
for parcelnum = 1:length(parcelIDs)
    parcel_withlines = surf_withlines .* (parcels_upsampled==parcelIDs(parcelnum));
    lines_within_parcel = unique(parcel_withlines);
    lines_within_parcel(lines_within_parcel==0)=[];
    
    parcel_communities = alt_consensusmap(parcelnum,:); parcel_communities(parcel_communities==0) = [];
    
    for comnum = 1:length(parcel_communities)
        
        lineval_indices = [comnum : length(parcel_communities) : length(lines_within_parcel)];
        linevals_for_this_community = lines_within_parcel(lineval_indices);
        for lineval = linevals_for_this_community'
            outputgifti(parcel_withlines==lineval) = parcel_communities(comnum);
        end
        
    end
    
end
save(gifti(single(outputgifti)),[outputfilestem '_confidence_164.func.gii']);
