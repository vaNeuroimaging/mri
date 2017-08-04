directory = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/Old/ICA_lowthresh20dim_HBMpaper.gica/ROIs/';

images = dir([directory '*.nii']);

for i = 1:length(images)
    o = maroi_image([directory images(i).name]);
    saveroi(o, [directory images(i).name(1:end-4) '_roi.mat']);
end