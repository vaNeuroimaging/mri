function r = spatial_correlation(Image1,Image2,Wholebrainmask)
%Spatial Correlation
% Conduct a spatial correlation between two MRI images to identify the
% degree of spatial similarity
%
%R = spatial_correlation('Image1','Image2','Wholebrainmask')
%
%R is the Pearson's R value returned by the spatial correlation.
%Image1 and Image2 are strings indicating the full-path locations of the
% two images to be compared.  These files should be in .hdr or .nii format.
%Wholebrainmask is a string indicating the full-path location of a
% binary whole-brain ROI within which the correlation will be tested (so
% that the correlation is not affected by irrelevant voxels outside the brain).
% Wholebrainmask may be a Marsbar ROI (in _roi.mat format), or it may be a
% .nii, .hdr, or .nii.gz file; if one of the latter, a _roi.mat file will be written
% out.
%
%Requires Marsbar to be in the Matlab path.
%
%Created by E. Gordon 11/04/09

filenamelength = max(length(Image1),length(Image2));

P = [sprintf(['%-' num2str(filenamelength) 's'],Image1); sprintf(['%-' num2str(filenamelength) 's'],Image2)];


if ~strcmp(Wholebrainmask(end-7:end),'_roi.mat')
    if strcmp(Wholebrainmask(end-6:end),'.nii.gz')
        FSL_from_matlab(['fslchfiletype NIFTI ' Wholebrainmask]);
        Wholebrainmask = Wholebrainmask(1:end-3);
    end
    
    
    outputname = [Wholebrainmask(1:end-3),'_roi.mat'];
    o = [];
    d = [];
    imgname = Wholebrainmask;
    [p f e] = fileparts(imgname);
    binf = 1;
    func = 'img >= 1';
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



brainroi = maroi(Wholebrainmask);

Y = getdata(brainroi,P,'l');

correlations = corrcoef(Y(1,:),Y(2,:));
r = correlations(2,1);


