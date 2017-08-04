tracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz');
tractvals = unique(tracts.img); tractvals(tractvals==0) = [];

tracts_to_noterode = [];%[11:14 21 22 35:40];

outtracts = tracts; outtracts.img(:) = 0;
temptracts = tracts;

for i = 1:length(tractvals)
    disp(num2str(i))
    temptracts.img = tracts.img .* uint8(tracts.img==tractvals(i));
    if ~any(tracts_to_noterode==tractvals(i))
        save_untouch_nii_2D(temptracts,'Temp.nii.gz')
        system('fslmaths Temp.nii.gz -kernel 3D -ero Temp_ero.nii.gz');
        temptracts = load_untouch_nii_2D('Temp_ero.nii.gz');
    end
    outtracts.img = outtracts.img + uint8(temptracts.img);
end

save_untouch_nii_2D(outtracts,'JHU-ICBM-labels-2mm-ero.nii.gz')