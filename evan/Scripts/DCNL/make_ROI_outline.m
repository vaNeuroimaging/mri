seedimagename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/TNN-aDMN-vmPFC_12mm_4space.nii';
seed = load_nii(seedimagename);
data = seed.img;
outdata = zeros(size(data));
for x = 2:size(data,1)-1
    for y = 2:size(data,2)-1
        for z = 2:size(data,3)-1
            if data(x,y,z) == 1 && mean(mean(mean(data(x-1:x+1,y-1:y+1,z-1:z+1))))<1
                outdata(x,y,z) = 1;
            end
        end
    end
end

seed.img = uint8(outdata);
save_nii(seed,'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/TNN-aDMN-vmPFC_12mm_4space_outline.nii');
        