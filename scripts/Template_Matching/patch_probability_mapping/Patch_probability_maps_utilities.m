%% Select maps of interest


ncortverts = 59412;

inds = [9 10 13 15:18 21 27 28 30 59 60 65 68 73:84 86 87 124:126 129 132 133 137 138 182 184 185 187 188 192:195 197 198 229 233 236 237 240 243 298:300 302:308 310 311 314 367 371 378 380 396 411 413 415 416 418 421:423 425 426 428 471 472 474 477 483 485 486 490 492 537 541 545 546 593 595 598:600 603 608 610 613 616 692 693 694 696 699 701 703 706 708 710 761 762:774 826:829 893 899 904 906 907 910];

mapname = 'Cluster_probability_maps_sorted_10mm';
sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines';
interestingfeature = 'notintemplate_plusextended';

removegroupcluster_inds  = [396;8];
dilate = 14;

groupsystems = ft_read_cifti_mod('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
%groupsystems = ft_read_cifti_mod('/data/cn4/evan/fsaverage_LR32k/Yeo_17_powercolors.dtseries.nii');
groupsystems = groupsystems.data(1:ncortverts);

baddata = ft_read_cifti_mod('/data/cn4/evan/ROIs/Baddata_bigcluster_LR_fillholes.dtseries.nii');
baddata = logical(baddata.data(1:ncortverts));

neighbors = cifti_neighbors('/data/cn4/evan/ROIs/Baddata_bigcluster_LR_fillholes.dtseries.nii');

distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat');
distances = distances(1:ncortverts,1:ncortverts);
radius = 5;



maps = ft_read_cifti_mod([mapname '.dscalar.nii']);
sizes = ft_read_cifti_mod([sizesname '.dscalar.nii']);
origsize = size(maps.data,1);
maps.data = maps.data(1:ncortverts,inds);
sizes.data = sizes.data(1:ncortverts,inds);
maps.mapname = maps.mapname(inds);
sizes.mapname = sizes.mapname(inds);



for i = 1:size(maps.data,2)
    
tokens = tokenize(maps.mapname{i},' ');
if strcmp(tokens{2}(end),':')
    tokens(2) = [];
end
thismap = maps.data(:,i) * str2num(tokens{6}) / str2num(tokens{4});
thissizes = sizes.data(:,i);

removenum = find(removegroupcluster_inds(1,:)==inds(i));

if ~isempty(removenum)
    ID = removegroupcluster_inds(2,removenum);

    removeinds = find(groupsystems==ID);
    dilate_distances = (distances(:,removeinds) < dilate);
    removeinds_dilated = logical(any(dilate_distances,2));
    
    dilate_distances = (distances(:,removeinds) < (dilate+5));
    removeinds_dilated_bigger = logical(any(dilate_distances,2));

    
    thismap(removeinds_dilated) = 0;
    thismap(baddata) = 0;
    
    tempmap = thismap; tempmap(removeinds_dilated_bigger) = 0;
    
    [ign,maxi] = max(tempmap);
    
    thismap(distances(:,maxi) > 80) = 0;
    
    temp = distances(:,maxi) < radius;
    tempneighs = unique(neighbors(logical(temp),2:end));
    tempneighs(isnan(tempneighs)) = [];
    borderinds = setdiff(tempneighs,find(temp));
    thissizes = zeros(size(thismap));
    thissizes(borderinds) = 1;
    
    maps.data(:,i) = thismap;
    sizes.data(:,i) = thissizes;
end
    
    



end

maps.data(end+1:origsize,:) = 0;
sizes.data(end+1:origsize,:) = 0;

ft_write_cifti_mod([mapname '_' interestingfeature '.dscalar.nii'],maps)
ft_write_cifti_mod([sizesname '_' interestingfeature '.dscalar.nii'],sizes)


%% Combine separated patch probability maps into one map per system

% mapname = 'Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108';
% sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108';

% mapname = 'Cluster_probability_maps_sorted_10mm_modified_108matchedto120';
% sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_modified_108matchedto120';

% mapname = 'Cluster_probability_maps_sorted_10mm_modified_HCPmatchedto120_matchedto108';
% sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_modified_HCPmatchedto120_matchedto108';

% mapname = 'Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108_renamed_matchedtoYeo_tweaked';
% sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108_renamed_matchedtoYeo_tweaked';

% mapname = 'Cluster_probability_maps_sorted_10mm_notintemplate_plusextended_renamed_YeomatchedtoPower_tweaked';
% sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_notintemplate_plusextended_renamed_YeomatchedtoPower_tweaked';

% mapname = 'Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108_inYeoTemplates';
% sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108_inYeoTemplates';

mapname = 'Cluster_probability_maps_sorted_10mm_notintemplate_plusextended_inPowerTemplates';
sizesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_notintemplate_plusextended_inPowerTemplates';

maps = ft_read_cifti_mod([mapname '.dscalar.nii']);
sizes = ft_read_cifti_mod([sizesname '.dscalar.nii']);

mapsout = maps;
sizesout = maps;
mapsout.data = zeros(size(maps.data,1),1);
mapsout.mapname = [];
sizesout.data = zeros(size(mapsout.data));
sizesout.mapname = [];

prevnetwork = [];
col = 0;

for i = 1:size(maps.data,2)
    
    if all(sum(maps.data(:,1:i-1) - repmat(maps.data(:,i),1,i-1),1))
    
    tokens = tokenize(maps.mapname{i},' ');
    newmap = maps.data(:,i) * str2num(tokens{6}) / str2num(tokens{4});
    newsizes = sizes.data(:,i);
    
    if ~strcmp(tokens{1},prevnetwork)
        col = col+1;
        prevnetwork = tokens{1};
        mapsout.data(:,col) = 0;
        sizesout.data(:,col) = 0;
        mapsout.mapname{col} = tokens{1}(1:end-1);
        sizesout.mapname{col} = tokens{1}(1:end-1);
    end
    
    mapsout.data(:,col) = mapsout.data(:,col) + newmap;
    sizesout.data(:,col) = sizesout.data(:,col) + newsizes;
    end
end

ft_write_cifti_mod([mapname '_combined.dscalar.nii'],mapsout)
ft_write_cifti_mod([sizesname '_combined.dscalar.nii'],sizesout)







%% Match to other datasets

ncortverts = 59412;

% mapname = 'Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended';
% outlinename = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended';
% basename = '120';
% 
% testmapsname = '/data/cn4/evan/RestingState/Ind_variability/HCP_allsys_clean30mm/Cluster_probability_maps_sorted_10mm_modified';
% testoutlinesname = '/data/cn4/evan/RestingState/Ind_variability/HCP_allsys_clean30mm/Cluster_probability_maps_sorted_10mm_mediansizeoutlines_modified';
% testname = 'HCP';

% mapname = 'Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended_matchedtoHCP';
% outlinename = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended_matchedtoHCP';
% basename = '120';
% 
% testmapsname = '/data/cn4/evan/RestingState/Ind_variability/108_allsys_clean30mm/Cluster_probability_maps_sorted_10mm_modified';
% testoutlinesname = '/data/cn4/evan/RestingState/Ind_variability/108_allsys_clean30mm/Cluster_probability_maps_sorted_10mm_mediansizeoutlines_modified';
% testname = '108';
% 
% bringalongmapsname = 'Cluster_probability_maps_sorted_10mm_modified_HCPmatchedto120';
% bringalongoutlinesname = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_modified_HCPmatchedto120';

mapname = 'Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108_renamed';
outlinename = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108_renamed';
basename = 'Power';

testmapsname = '/data/cn4/evan/RestingState/Ind_variability/120_108_YeoTemplate_betterkden_clean30mm/120/Cluster_probability_maps_sorted_10mm_notintemplate_plusextended_renamed';
testoutlinesname = '/data/cn4/evan/RestingState/Ind_variability/120_108_YeoTemplate_betterkden_clean30mm/120/Cluster_probability_maps_sorted_10mm_mediansizeoutlines_notintemplate_plusextended_renamed';
testname = 'Yeo';



maps = ft_read_cifti_mod([mapname '.dscalar.nii']); 
outlines = ft_read_cifti_mod([outlinename '.dscalar.nii']);

testmaps = ft_read_cifti_mod([testmapsname '.dscalar.nii']);
testmaps.data = testmaps.data(1:ncortverts,:); 
testmaps.data((ncortverts+1):size(maps.data,1),:) = 0;

testoutlines = ft_read_cifti_mod([testoutlinesname '.dscalar.nii']);
testoutlines.data = testoutlines.data(1:ncortverts,:); 
testoutlines.data((ncortverts+1):size(maps.data,1),:) = 0;


tokens = tokenize(testmapsname,'/');
testmapswritename = tokens{end};
tokens = tokenize(testoutlinesname,'/');
testoutlineswritename = tokens{end};


thresh = .05;

mapsout = maps;
mapsout.data = [];
mapsout.mapname = [];
outlinesout = mapsout;
matchedmaps = mapsout;
matchedoutlines = mapsout;

if exist('bringalongmapsname')
    bringalongmaps = ft_read_cifti_mod([bringalongmapsname '.dscalar.nii']);
    bringalongoutlines = ft_read_cifti_mod([bringalongoutlinesname '.dscalar.nii']);
    
    bringalongmapsout = bringalongmaps;
    bringalongmapsout.data = [];
    bringalongmapsout.mapname = [];
    bringalongoutlinesout = bringalongoutlines;
    bringalongoutlinesout.data = [];
    bringalongoutlinesout.mapname = [];
    
    tokens = tokenize(bringalongmapsname,'/');
    bringalongmapswritename = tokens{end};
    tokens = tokenize(bringalongoutlinesname,'/');
    bringalongoutlineswritename = tokens{end};
end


col = 0;

for i = 1:size(maps.data,2)
    tokens = tokenize(maps.mapname{i},' ');
    mapnetwork = tokens{1};
    mappct = str2num(tokens{6}) / str2num(tokens{4});
    correls = zeros(1,size(testmaps.data,2));
    for j = 1:size(testmaps.data,2);
        tokens = tokenize(testmaps.mapname{j},' ');
        testmapnetwork = tokens{1};
        testmappct = str2num(tokens{6}) / str2num(tokens{4});
        if strcmp(mapnetwork,testmapnetwork) %&& (abs(mappct - testmappct) < pctdiffthresh)
            inds = (maps.data(:,i)>=.05) | (testmaps.data(:,j)>=.05);
            correls(j) = paircorr_mod(maps.data(inds,i),testmaps.data(inds,j));
        end
    end
    [maxcorrel(i), maxi] = max(correls);
    if maxcorrel(i) > thresh
        col = col+1;
        mapsout.data(:,col) = maps.data(:,i);
        mapsout.mapname{col} = maps.mapname{i};
        outlinesout.data(:,col) = outlines.data(:,i);
        outlinesout.mapname{col} = outlines.mapname{i};
        matchedmaps.data(:,col) = testmaps.data(:,maxi);
        matchedmaps.mapname{col} = testmaps.mapname{maxi};
        matchedoutlines.data(:,col) = testoutlines.data(:,maxi);
        matchedoutlines.mapname{col} = testoutlines.mapname{maxi};
        
        if exist('bringalongmapsname')
            bringalongmapsout.data(:,col) = bringalongmaps.data(:,i);
            bringalongmapsout.mapname{col} = bringalongmaps.mapname{i};
            bringalongoutlinesout.data(:,col) = bringalongoutlines.data(:,i);
            bringalongoutlinesout.mapname{col} = bringalongoutlines.mapname{i};
        end
    end
end
        
ft_write_cifti_mod([mapname '_matchedto' testname],mapsout)
ft_write_cifti_mod([outlinename '_matchedto' testname],outlinesout)
ft_write_cifti_mod([testmapswritename '_' testname 'matchedto' basename],matchedmaps)
ft_write_cifti_mod([testoutlineswritename '_' testname 'matchedto' basename],matchedoutlines)

if exist('bringalongmapsname')
    ft_write_cifti_mod([bringalongmapsname '_matchedto' testname],bringalongmapsout)
    ft_write_cifti_mod([bringalongoutlinesname '_matchedto' testname],bringalongoutlinesout)
end
    



%% Regions from median size outlines


outlinesfile = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108';
neighbors = cifti_neighbors([outlinesfile '.dscalar.nii']);
outlines = ft_read_cifti_mod([outlinesfile '.dscalar.nii']);

dtseries_template = ft_read_cifti_mod('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');

distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat');

networknames = {'Default','Visual','FrontoPar','DorsalAttn','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'};
networknums = [1 2 3 5 7 8 9 10 11 12 13 14 15 16];

outlines_combined = zeros(size(outlines.data,1),1);
for i = 1:size(outlines.data,2);
    tokens = tokenize(outlines.mapname{i},' ');
    networkind = strcmp(tokens{1}(1:end-1),networknames);
    networknum = networknums(networkind);
    
    thisoutline = outlines.data(:,i);
    thisoutline_inverted = ~thisoutline;
    inds_toremove = all(distances(:,logical(thisoutline))>50,2);
    inds_toremove(59413:end) = 1;
    thisoutline_inverted(inds_toremove) = 0;
    clusters = cifti_cluster(thisoutline_inverted,.5,1.5,0,neighbors);
    
    [~,sorti] = sort(sum(clusters,1),'descend');
    secondbiggest = sorti(2);
    
    this_cluster = clusters(:,secondbiggest) * networknum;
    
    
    overlap_withprev = logical(outlines_combined) & logical(this_cluster);
    outlines_combined = outlines_combined + this_cluster;
    outlines_combined(overlap_withprev) = -1;

end

dtseries_template.data = outlines_combined;
ft_write_cifti_mod([outlinesfile '_all'],dtseries_template)


%% Size, frequency, and spatial correlations of matched patches

patches = {'Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108','Cluster_probability_maps_sorted_10mm_modified_108matchedto120','Cluster_probability_maps_sorted_10mm_modified_HCPmatchedto120_matchedto108'};
names = {'Wash U','Dartmouth','HCP'};

ignoreinds = [1 3 13 14 30 42 65];




for i = 1:length(patches)
    patchesi = ft_read_cifti_mod([patches{i} '.dscalar.nii']);
    
    if i==1
        
        maps = [1:size(patchesi.data,2)];
        maps = setdiff(maps,ignoreinds);
        
        sizes = zeros(length(maps),length(patches));
        subpcts = zeros(length(maps),length(patches));
    end
    
    for mapnum = 1:length(maps);
        tokensi = tokenize(patchesi.mapname{maps(mapnum)},' ');
        subpcts(mapnum,i) = str2num(tokensi{6}) / str2num(tokensi{4});
        sizes(mapnum,i) = str2num(tokensi{10}(6:end-4));
    end
    
    
    
    
    for j = (i+1):length(patches)
        patchesj = ft_read_cifti_mod([patches{j} '.dscalar.nii']);
        
        spatialcorrels{i,j} = zeros(length(maps),1);
        
        for mapnum = 1:length(maps);
            tokensi = tokenize(patchesi.mapname{maps(mapnum)},' ');
            
            tokensj = tokenize(patchesj.mapname{maps(mapnum)},' ');
            
            patchi_allpct = patchesi.data(:,maps(mapnum)) * str2num(tokensi{6}) / str2num(tokensi{4});
            patchj_allpct = patchesj.data(:,maps(mapnum)) * str2num(tokensj{6}) / str2num(tokensj{4});
            inds = (patchi_allpct > .05) | (patchj_allpct > .05);
            spatialcorrels{i,j}(mapnum) = paircorr_mod(patchi_allpct(inds),patchj_allpct(inds));
        end
    end
    
end


    
%% Table of patch locations

coordsfile = '/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_MNIcoords.dtseries.nii';
MNIcoords = ft_read_cifti_mod(coordsfile); MNIcoords = MNIcoords.data;

labelnums = ft_read_cifti_mod('/data/cn4/evan/ROIs/mode_LR_aparc_a2009s_32k_fs_LR.dtseries.nii');
labelnums = labelnums.data;

labels = gifti('/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc33416/7112b_fs_LR/fsaverage_LR32k/vc33416.L.aparc.a2009s.32k_fs_LR.label.gii');
labelnames =  labels.private.label.name;
labelkey = labels.private.label.key;

maps = ft_read_cifti_mod('Cluster_probability_maps_sorted_10mm_40sub_notintemplate_plusextended.dscalar.nii');


outputfilename = 'Probability_maps_notintemplate_table.txt';
delete([outputfilename]);
fid = fopen([outputfilename],'at'); %open the output file for writing
fclose(fid);


for mapnum = 1:size(maps.data,2)
    
    [~,maxind] = max(maps.data(:,mapnum));
    
    if maxind <= 29696
        hem = 'left';
    else
        hem = 'right';
    end
    
    maxcoords = [num2str(round(MNIcoords(maxind,1))) ' ' num2str(round(MNIcoords(maxind,2))) ' ' num2str(round(MNIcoords(maxind,3)))];
    
    maxlabelnum = labelnums(maxind);
    maxlabel = labelnames{labelkey==maxlabelnum};
    
    tokens = tokenize(maps.mapname{mapnum},' ');
    
    system = tokens{1}(1:end-1);
    nsubs = [num2str(round(str2num(tokens{6}) / str2num(tokens{4}) * 100)) '%'];
    mapsize = tokens{10}(6:end-1);
    
    texttowrite = sprintf('%s\t%s\t%s\t%s\t%s\t%s',system,nsubs,mapsize,hem,maxcoords,maxlabel);
    
    dlmwrite(outputfilename,texttowrite,'-append','delimiter','')
end


