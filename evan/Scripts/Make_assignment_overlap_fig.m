ParcelID = 16715;
outputfilename = ['/data/cn4/evan/RestingState/Ind_variability/Subjects/Target230/Assignment_overlap_parcel' num2str(ParcelID) '.func.gii'];

Parcelfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_watershedmerge.func.gii';
parcels = gifti(Parcelfilename); parcels = parcels.cdata;

Groupmatchfilename = '/data/cn4/evan/RestingState/Ind_variability/Subjects/Target230/Group_parcels_assignedtosubs_L.func.gii';
groupmatch = gifti(Groupmatchfilename); groupmatch = groupmatch.cdata;

Submatchfilename = '/data/cn4/evan/RestingState/Ind_variability/Subjects/Target230/Subject_parcels_assignedtogroup_L.func.gii';
submatch = gifti(Submatchfilename); submatch = submatch.cdata;


bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


outputfile = zeros(size(parcels));

%outputfile(parcels==ParcelID) = 1;

for s = 1:size(submatch,2)
    
    matchID = mean(groupmatch(parcels==ParcelID,s));
    
    subparcelindices = find(submatch(:,s)==matchID);
    
    outlineindices = [];
    
    for index = subparcelindices'
        
        indneighs = neighbors(index,2:7);
        indneighs(isnan(indneighs))=[];
        
        if any(submatch(indneighs,s)==0)
            
            outlineindices = [outlineindices index];
        end
    end
    
%     outlinetracker = ones(1,length(outlineindices));
%     outlinetracker(1) = 0;
%     
%     for indexnum = 2:length(outlineindices)
%         index = outlineindices(indexnum);
%         
%         indneighs = neighbors(index,2:7);
%         indneighs(isnan(indneighs))=[];
%         outlineneighs = intersect(outlineindices(logical(outlinetracker)),indneighs);
%         
%         if length(outlineneighs) > 1;
%             outlinetracker(indexnum) = 0;
%         end
%     end
%     
%     outlineindices(outlinetracker==0) = [];
    
    outputfile(outlineindices) = s+1;
            
    
end

save(gifti(single(outputfile)),outputfilename);
    

