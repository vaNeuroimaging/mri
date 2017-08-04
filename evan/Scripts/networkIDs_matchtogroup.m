function networkIDs_matchtogroup(subnetworksfilename,distances)
%networkIDs_matchtogroup(subnetworksfilename,[distances])

slashlocs = strfind(subnetworksfilename,'/');
if isempty(slashlocs)
    subnetworksfolder = pwd;
    subnetworksfile = subnetworksfilename;
else
    subnetworksfolder = subnetworksfilename(1:slashlocs(end));
    subnetworksfile = subnetworksfilename(slashlocs(end)+1 : end);
end


groupnetworksfile = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR.dtseries.nii';

nsurfverts = 29696 + 29716;


if ~exist('distances')
    load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
    distances = distances(1:nsurfverts,1:nsurfverts);
end

extraassignments = [11.5 4 6 17 2.5 4.5 18:100]; 

subnetworks = cifti_read(subnetworksfile);
surfsubnetworks = subnetworks(1:nsurfverts);
subIDs = unique(surfsubnetworks); subIDs(subIDs<1) = [];

groupnetworks = cifti_read(groupnetworksfile);
surfgroupnetworks = groupnetworks(1:nsurfverts);
groupIDs = unique(surfgroupnetworks); groupIDs(groupIDs<1) = [];

matched = zeros(size(surfsubnetworks));


mean_mindists = zeros(length(groupIDs),length(subIDs));
for groupIDnum = 1:length(groupIDs)
    groupID = groupIDs(groupIDnum);
    %lowest_mindist = Inf;
    %mindist_ID = 0;

    
    for subIDnum = 1:length(subIDs)
        subID = subIDs(subIDnum);
        submindists = min(distances(surfgroupnetworks==groupID,surfsubnetworks==subID),[],1);
        submindists(submindists > 50) = [];
        groupmindists = min(distances(surfgroupnetworks==groupID,surfsubnetworks==subID),[],2);
        groupmindists(groupmindists>50) = [];
        mean_mindists(groupIDnum,subIDnum) = mean([submindists(:); groupmindists(:)]);
        
    %    if mean_mindists < lowest_mindist
    %        lowest_mindist = mean_mindists;
    %        mindist_ID = subID;
        %end
    end
    %unassignedIDs(unassignedIDs==mindist_ID) = [];
    %matched(surfsubnetworks==mindist_ID) = groupID;
end

unassignedIDs = subIDs;
[assign,ign] = munkres(mean_mindists);
for assignednum = 1:length(assign)
    if (assign(assignednum)) > 0 && (mean_mindists(assignednum,assign(assignednum)) < 10)
        matched(surfsubnetworks==subIDs(assign(assignednum))) = groupIDs(assignednum);
        unassignedIDs(unassignedIDs==subIDs(assign(assignednum))) = [];
    end
end

for extraIDnum = 1:length(unassignedIDs)
    matched(surfsubnetworks==unassignedIDs(extraIDnum)) = extraassignments(extraIDnum);
end

out = zeros(size(groupnetworks));
out(1:nsurfverts) = matched;
outfilename = [subnetworksfolder '/Powercolors_' subnetworksfile];
cifti_write_wHDR(out,[],outfilename)