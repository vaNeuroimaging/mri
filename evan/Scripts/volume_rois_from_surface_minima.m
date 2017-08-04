minimafile = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_minima.func.gii';
surfmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT.func.gii';
labelsfile = '/data/cn4/evan/ROIs/FinalLabels/Power_Neuron11_dil.L.32k_fs_LR.label.gii';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/1mmROIs';
ROIradius = 1;
coordsfile = '/data/cn4/evan/fsaverage_LR32k/node_coords.txt';
maskfile='/data/cn4/evan/ROIs/glm_atlas_mask_333.nii';
volcenter = [24.5 36 21];
volvoxsize = 3;


minimadata = gifti(minimafile);
minimadata = minimadata.cdata;

if ~isempty(surfmaskfile)

maskdata = gifti(surfmaskfile);
maskdata = maskdata.cdata;

minimadata = minimadata .* maskdata;

end

labelsdata = gifti(labelsfile);
labelsdata = labelsdata.cdata;

nodecoords = dlmread(coordsfile);

networkcounter = [];

for minimaloc = find(minimadata)'

    %dlmwrite([outputfolder 'minimalocations.txt'],num2str(minimaloc),'-append','delimiter','');
    
    networknum = labelsdata(minimaloc);
    
    if length(networkcounter) < networknum
        networkcounter(networknum) = 0;
    end
    networkcounter(networknum) = networkcounter(networknum) +1;
    
    [ign networkstring] = system(['awk ''/<Label Key="' num2str(networknum) '"/{print}'' ' labelsfile]);
    stringindex1 = max([strfind(networkstring,['<![CDATA[u' num2str(networknum)]) strfind(networkstring,['<![CDATA[a' num2str(networknum)])]) + length(['<![CDATA[u' num2str(networknum) '_']);
    stringindex2 = strfind(networkstring,']]') - 1;
    networkname = [networkstring(stringindex1:stringindex2) '_' num2str(networkcounter(networknum))];
    
    fslcoords = nodecoords(minimaloc,1:3) / volvoxsize + volcenter -1;
    
    system(['fslmaths ' maskfile ' -roi ' num2str(round(fslcoords(1))) ' 1 ' num2str(round(fslcoords(2))) ' 1 ' num2str(round(fslcoords(3))) ' 1 0 1 ' outputfolder filesep 'minimum_' networkname ' -odt float']);
    system(['fslmaths ' outputfolder filesep 'minimum_' networkname ' -kernel sphere ' num2str(ROIradius) ' -fmean -bin -mul ' maskfile ' ' outputfolder filesep 'minimum_' networkname ' -odt float']);
    system(['fslchfiletype NIFTI ' outputfolder filesep 'minimum_' networkname '.nii.gz']);

end
    
    
    
    