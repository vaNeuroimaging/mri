assignmentsfile = ['rawassn.txt'];
parcelfilename = {'/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_L_watershedmerge_0.6_tweaked.func.gii','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_R_watershedmerge_0.6_tweaked.func.gii'};
%{'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_L_nosmooth_watershedmerge_0.4_tweaked.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_R_nosmooth_watershedmerge_0.4_tweaked.func.gii'}; %{'/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_edge_L_watershedmerge_0.45_tweaked.func.gii','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_subsurf_edge_R_watershedmerge_0.45_tweaked.func.gii'};
templatefile = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';
outputfilestem = 'Poldrome_subsurf_LR_minsize5_consensus';
minsize = 5;

mincol = 10;

recolor = [];
%recolor = [1 10; 2 9; 3 3; 4 1; 5 7; 6 4.5; 7 6; 8 16; 9 5; 10 15; 11 2; 12 2.5; 13 11.5; 14 12;15 17];

power_surf_colormap = [.67 .67 .67;.67 .67 .67;1 0 0;0 0 .6;1 1 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;1 1 .8;0 .4 0;.25 .25 .25];


%% Size threshold

simplified = modify_clrfile('simplify',assignmentsfile,minsize);


%% Regularize

regularized = rawoutput2clr(simplified);
regularized(regularized<=1) = 0; regularized = regularized-1;
dlmwrite([assignmentsfile(1:end-4) '_minsize' num2str(minsize) '_regularized.txt'],regularized,'delimiter','\t');


if max(unique(regularized)) < (size(power_surf_colormap,1)-2)
    power_surf_colormap_touse = power_surf_colormap(1:max(unique(regularized))+2,:);
else
    power_surf_colormap_touse = [power_surf_colormap ; repmat(power_surf_colormap(end,:),(max(unique(regularized)) - (size(power_surf_colormap,1)-2)),1)];
end


%% Create initial consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments

consensusmap = regularized(:,mincol);


unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = regularized(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end


%% Reapply size threshold

consensusvals = unique(consensusmap);
for val = consensusvals(:)'
    if nnz(consensusmap==val) < minsize
        toosmallinds = find(consensusmap==val);
        for smallindex = toosmallinds(:)'
            thisassignments = regularized(unassignedindex,mincol:end);
            thisassignments(thisassignments<1) = [];
            thisassignments(thisassignments==val) = [];
            if ~isempty(thisassignments)
                consensusmap(smallindex) = thisassignments(1);
            end
        end
        
    end
end

% %% Move values around to make things prettier
% 
% vals_inconsensus = unique(consensusmap);
% morevals_in_matrix = setdiff(unique(regularized),vals_inconsensus);
% for val = 1:length(morevals_in_matrix)
%     regularized(regularized==morevals_in_matrix(val)) = max(vals_inconsensus) + val;
% end
% 
% for val = 1:max(vals_inconsensus);
%     if (nnz(consensusmap==val)==0) && (val<max(vals_inconsensus))
%         consensusmap(consensusmap==max(vals_inconsensus)) = val;
%         regularized(regularized==max(vals_inconsensus)) = val;
%         vals_inconsensus = unique(consensusmap);
%     end
% end




%% Recolor

temp_consensusmap = consensusmap;
temp_regularized = regularized;
for i = 1:size(recolor,1)
    temp_consensusmap(consensusmap==recolor(i,1)) = recolor(i,2);
    temp_regularized(regularized==recolor(i,1)) = recolor(i,2);
end
consensusmap = temp_consensusmap;
regularized = temp_regularized;


%% Write parcel output

figure
imagesc(sortrows(regularized,[size(regularized,2):-1:1]))
colormap(power_surf_colormap_touse)
title('Assignments')


templatedata = gifti(templatefile);


outputdata = zeros(size(templatedata.cdata,1),1);
crossthresh_outputdata = zeros(size(templatedata.cdata,1),size(regularized,2));

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
            crossthresh_outputdata(find(parcels==parcelIDs(parcelnum)) + (strcmp('R',hem)*ncortexLverts),col) = regularized(parcelnum + (strcmp('R',hem)*nLparcelIDs),col);
        end
    end
end
dlmwrite([outputfilestem '.txt'],consensusmap)
cifti_write_wHDR(outputdata,templatefile,outputfilestem)
cifti_write_wHDR(crossthresh_outputdata,templatefile,[outputfilestem '_allthresh'])
