hems = {'L','R'};
cd /data/cn4/evan/ROIs/Craddock/
for i = 4:42
    gunzip(['tcorr05_2level_all' sprintf('%04i',i) '.nii.gz'])
    map_vol_to_surface(['tcorr05_2level_all' sprintf('%04i',i) '.nii'],'both','enclosing')
    
    for hemnum = 1:2
        hem = hems{hemnum};
        
        Fill_Craddock_parcel_holes(['tcorr05_2level_all' sprintf('%04i',i) '_' hem '.func.gii'],['Craddock' sprintf('%04i',i) '_' hem '.func.gii'],hem)
        load(['/data/cn4/evan/RestingState/FC_Mapping_120/cov_corr_normalwall_' hem '.mat'])
        eval(['cov_corr = cov_corr_' hem '; clear cov_corr_' hem])
        generate_rotated_parcels_andPCA6(['Craddock' sprintf('%04i',i) '_' hem '.func.gii'],100,cov_corr,1,hem,['Craddock' sprintf('%04i',i) '_' hem '_'])
        
    end
end