
ROIpath = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/Cluster/';


imagestowrite = dir([ROIpath '*.nii']);


for filenum = 1:length(imagestowrite)
    
    imgname = [ROIpath imagestowrite(filenum).name];
    
    outputname = [imgname(1:end-4) '_roi.mat'];
    
    o = [];
    d = [];
    
    [p f e] = fileparts(imgname);
    binf = 1;
    
    func = 'img > 0';
    
    d = f; l = f;
    if ~isempty(func)
        d = [d ' func: ' func];
        l = [l '_f_' func];
    end
    if binf
        d = [d ' - binarized'];
        l = [l '_bin'];
    end
    o = maroi_image(struct('vol', spm_vol(imgname), 'binarize',binf,...
        'func', func));
    
    % convert to matrix format to avoid delicacies of image format
    o = maroi_matrix(o);
    
    
    o = descrip(o,d);
    o = label(o,l);
    
    varargout = {saveroi(o, outputname)};
    
end

