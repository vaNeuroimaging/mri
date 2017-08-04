tempnames = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Power_ROIs/MNI position names.txt','%s','delimiter','\t');

for i = 3:length(tempnames)
    names(ceil(i/2)-1,(2-rem(i,2))) = tempnames(i);
end

networknames = {'Limb','CO','Cereb','FPCl','SM','pDMN','DA','Vis','aDMN','FT','FPCr','Aud'};
networknums = [1 3 4 6 7 8 10 12 13 14 16 19];

allnetworks = load_nii('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/Old/ICA_lowthresh20dim_HBMpaper.gica/groupmelodic.ica/melodic_IC.nii');

for network = 1:length(networknames)
    
    networkmat(:,:,:,network) = allnetworks.img(:,:,:,networknums(network));
    
end

for ROInum = 1:size(names,1);
    xs = findstr('x',names{ROInum,1});
    negmnicoords = [-str2num(names{ROInum,1}(1:xs(1)-1)) str2num(names{ROInum,1}(xs(1)+1:xs(2)-1)) str2num(names{ROInum,1}(xs(2)+1:end))];
    matrixcoords = round(negmnicoords / [allnetworks.hdr.hist.srow_x(1:3);allnetworks.hdr.hist.srow_y(1:3);allnetworks.hdr.hist.srow_z(1:3)] + allnetworks.hdr.hist.originator(1:3));
    [maxvalue, networkindex] = max(networkmat(matrixcoords(1),matrixcoords(2),matrixcoords(3),:));
    disp([names{ROInum,1} ', ' names{ROInum,2} ': ' num2str(maxvalue) ' in ' networknames{networkindex}])
    
    origname = [names{ROInum,1}(1:xs(1)-1) '_' names{ROInum,1}(xs(1)+1:xs(2)-1) '_' names{ROInum,1}(xs(2)+1:end) '_roi.mat'];
    newname = [networknames{networkindex} '_' names{ROInum,2} '_roi.mat'];
    newname(strfind(newname,' ')) = '_';
    
    try
        copyfile(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Power_ROIs/Original/' origname],['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Power_ROIs/' newname]);
    catch
        disp('Failed')
    end
end