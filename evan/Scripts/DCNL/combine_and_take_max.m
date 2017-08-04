basedir = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/';

dirstocombine = {'Negative_DL1_Rest','Negative_DL2_Rest','Negative_DL3_Rest','Negative_DL4_Rest','Negative_DL5_Rest','Negative_DL6_Rest'};

for imagenum = 1:length(dirstocombine)
    thisimage = load_nii([basedir dirstocombine{imagenum} '/spmT_0001.hdr']);
    bigimage(:,:,:,imagenum) = thisimage.img;
end

maximage = max(bigimage,[],4);

for imagenum = 1:length(dirstocombine)
    imagetowrite = thisimage;    
    imagetowrite.img = squeeze(bigimage(:,:,:,imagenum)) .* (squeeze(bigimage(:,:,:,imagenum))==maximage);
    imagetowrite.fileprefix = [basedir dirstocombine{imagenum} '/spmT_0001_max.hdr'];
    save_nii(imagetowrite,[basedir dirstocombine{imagenum} '/spmT_0001_max.hdr']);
    clear imagetowrite
end