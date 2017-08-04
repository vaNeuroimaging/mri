 parcel_filename = ['/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_R_crossthresh_watershedmerge.func.gii'];
%parcel_filename = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_L_crossthresh_watershedmerge_tweaked.func.gii';
outputfile = '/data/cn4/evan/RestingState/Consensus_parcels/Poldrome_R_Final.func.gii';
assignments = Final;

recolor = [];
%recolor = [1 2;2 10;3 9;4 7;5 12;6 11;7 16;8 3;9 5;10 1;11 4];
%recolor = [1 2; 2 10;3 5;4 9;5 1;6 7;7 21;8 13;9 12;10 3;11 20;12 18;13 17;14 14;15 11;16 21];
temp_assignments = assignments;
for i = 1:size(recolor,1)
    temp_assignments(assignments==recolor(i,1)) = recolor(i,2);
end
assignments = temp_assignments;

parcels = gifti(parcel_filename); parcels = parcels.cdata;
parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0) = [];

outputdata = zeros(size(parcels,1),size(assignments,2));

for colnum = 1:size(assignments,2)
    for parcelnum = 1:length(parcelIDs)
        
        outputdata(parcels==parcelIDs(parcelnum),colnum) = assignments(parcelnum,colnum);
        
    end
end

save(gifti(single(outputdata)),outputfile)