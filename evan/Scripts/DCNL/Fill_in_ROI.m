imgname = '/fmri/data3/Evan/Gene-Rest-Nback/Data/101/Struct/FIRST_Hip_seg/101_VT.nii.gz';

if strcmp(imgname(end-2:end),'.gz')
    eval(['!fslchfiletype NIFTI ' imgname ])
    imgname = imgname(1:end-3);
end

image = load_nii(imgname);
image.img(find(image.img)) = 1;
newimage = image;

for y = 1:size(image.img,2);
    [xcoord, zcoord] = find(squeeze(image.img(:,y,:)));
    for i = 1:length(xcoord)
        newimage.img(xcoord(i),y-1,zcoord(i)) = 1;
        newimage.img(xcoord(i),y+1,zcoord(i)) = 1;
    end
end

save_nii(newimage,'/fmri/data3/Evan/Gene-Rest-Nback/Data/101/Struct/FIRST_Hip_seg/101_VT_filledin.nii');