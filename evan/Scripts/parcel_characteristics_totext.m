parcels = cifti_read('parcelnum.dtseries.nii');
assignments = cifti_read('120_108_combined_LR_infomap_in120/watershed_Tk0005to005in0001_S1to1_surfxd30_INFMAP/tallyconsensus_120_108_LR_recolored2.dtseries.nii');
IDs = unique(parcels); IDs(IDs==0) = [];

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
ncortvertsL = nnz(maskL);

sizesL = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR_surfaceareas.func.gii');
sizesL = sizesL.cdata(logical(maskL));
sizesR = gifti('/data/cn4/evan/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR_surfaceareas.func.gii');
sizesR = sizesR.cdata(logical(maskR));

outputfilename = '/data/cn5/selfRegulation/V4Process_nosmooth/gradients_120_108_combined_subsurf_nosmooth/Parcels.txt';

delete([outputfilename]);
fid = fopen([outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','ParcelID','Hem','Size','Community','Name'); %write the output file header
fclose(fid);
dlmwrite([outputfilename],' ','-append');

parcelcommunityfilename = '/data/cn5/selfRegulation/V4Process_nosmooth/gradients_120_108_combined_subsurf_nosmooth/ParcelCommunities.txt';
fid = fopen(parcelcommunityfilename,'at');
fclose(fid);

%parcelsout = zeros(length(parcels),length(IDs));


structL = gifti('/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii'); structL = structL.cdata(logical(maskL));
structR = gifti('/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii'); structR = structR.cdata(logical(maskR));

networknames = {'None','None','Default','Visual','FrontoParietal','FPC2','DorsalAttn','DorsAttn2','VentralAttn','Salience','CinguloOperc','SMhand','SMmouth','Auditory','MTL1','MTL2','PERN','Context'};

bufsize = 524288;
headertext = textread('/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc33416/7112b_fs_LR/fsaverage_LR32k/vc33416.L.aparc.a2009s.32k_fs_LR.label.gii','%s','delimiter','\r','bufsize',bufsize);

count = 1;

for row = 1:length(headertext)
    if length(headertext{row}) > 15 && strcmp(headertext{row}(end-10:end),']]></Label>')
        %if count > 0
            beginfind = strfind(headertext{row},'[CDATA[') + 7;
            
            thislabel = headertext{row}(beginfind : end-11);
            
            if length(thislabel) > 8 && strcmp(thislabel(1:8),'G_and_S_')
                 label{count} = thislabel(9:end);
            elseif strcmp(thislabel(1:2),'G_')
                label{count} = thislabel(3:end);
                %label{count} = [thislabel(3:end) '_gyrus'];
            elseif strcmp(thislabel(1:2),'S_')
                label{count} = thislabel(3:end);
                %label{count} = [thislabel(3:end) '_sulcus'];
            else
                label{count} = thislabel;
            end
            
        %end
        count = count+1;
    end
end



for IDnum = 1:length(IDs)
    inds = find(parcels==IDs(IDnum));

    %parcelsout(inds,IDnum) = 1;
    
    community = mean(assignments(inds));
    
    if all(inds<=ncortvertsL)
        hem = 'L';
        size = sum(sizesL(inds));
        modestruct = mode(structL(inds));
    else
        hem = 'R';
        size = sum(sizesR(inds - ncortvertsL));
        modestruct = mode(structR(inds - ncortvertsL));
    end
    
    count = 1;
    name{IDnum} = [networknames{community+2} '_' hem '_' label{modestruct+1} '_' num2str(count)];
    
    for i = 1:IDnum-1
        if strcmp(name{i},name{IDnum})
            count = count+1;
            name{IDnum} = [networknames{community+2} '_' hem '_' label{modestruct+1} '_' num2str(count)];
        end
    end
        
    
    texttowrite = [num2str(IDs(IDnum)) ' ' hem ' ' num2str(size) ' ' networknames{community+2} ' ' name{IDnum}];
    
    dlmwrite(outputfilename,texttowrite,'delimiter','','-append')
    
    dlmwrite(parcelcommunityfilename,community,'-append')
end

%cifti_write_wHDR(parcelsout,[],'parcelnum_splitup')
    
    