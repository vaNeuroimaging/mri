minimafolder = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/voxelwiseFC/';
surfmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT.func.gii';
%labelsfile = '/data/cn4/evan/ROIs/FinalLabels/Power_Neuron11_dil.L.32k_fs_LR.label.gii';
outputfolder = minimafolder;
ROIradius = 4;
coordsfile = '/data/cn4/evan/fsaverage_LR32k/node_coords.txt';
maskfile='/data/cn4/evan/ROIs/glm_atlas_mask_333.nii';
volcenter = [24.5 36 21];
volvoxsize = 3;


minimafiles = dir([minimafolder '/minima_Avg*.func.gii']);
for minimanum = 1:length(minimafiles)
    thisminima = gifti([minimafolder '/' minimafiles(minimanum).name]);
    allminimamatrix(:,minimanum) = thisminima.cdata;
end

minimadata = double(logical(sum(allminimamatrix,2)));


if ~isempty(surfmaskfile)

maskdata = gifti(surfmaskfile);
maskdata = maskdata.cdata;

minimadata = minimadata .* (~maskdata);

end


nodecoords = dlmread(coordsfile);

networkcounter = [];

for minimaloc = find(minimadata)'

    %dlmwrite([outputfolder 'minimalocations.txt'],num2str(minimaloc),'-append','delimiter','');
    
    orig_minimanum = find(allminimamatrix(minimaloc,:));
    orig_minimanum = orig_minimanum(1);
    
    if length(networkcounter) < orig_minimanum
        networkcounter(orig_minimanum) = 0;
    end
    networkcounter(orig_minimanum) = networkcounter(orig_minimanum) +1;

    networkname = [minimafiles(orig_minimanum).name(1:end-9) '_' num2str(networkcounter(orig_minimanum))];
    
    fslcoords = nodecoords(minimaloc,1:3) / volvoxsize + volcenter -1;
    
    system(['fslmaths ' maskfile ' -roi ' num2str(round(fslcoords(1))) ' 1 ' num2str(round(fslcoords(2))) ' 1 ' num2str(round(fslcoords(3))) ' 1 0 1 ' outputfolder filesep 'ROI_' networkname ' -odt float']);
    system(['fslmaths ' outputfolder filesep 'ROI_' networkname ' -kernel sphere ' num2str(ROIradius) ' -fmean -bin -mul ' maskfile ' ' outputfolder filesep 'ROI_' networkname ' -odt float']);
    system(['fslchfiletype NIFTI ' outputfolder filesep 'ROI_' networkname '.nii.gz']);

end
    
    
    
    