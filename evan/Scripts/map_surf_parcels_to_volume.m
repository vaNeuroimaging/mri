dir = '/data/cn5/selfRegulation/V4Process_nosmooth/gradients_120_108_combined_subsurf_nosmooth/volumeparcels/';
cd(dir)
HEMS = {'L';'R'};
surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
outtextname = 'ParcelCommunities.txt';
dlmwrite(outtextname,'')
for h = 1:2
    
        water = gifti([dir '/../parcelnum_' HEMS{h} '.func.gii']);
        water = water.cdata;
        
        %communities = gifti(['/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_LR_infomap/watershed_Tk0005to005in0001_S1to1_surfxd20_INFMAP/120_subsurf_' HEMS{h} '_minsize5_consensus.func.gii']);
        %communities = communities.cdata;

    mask = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{h} '.32k_fs_LR.func.gii']);
    mask = mask.cdata;
    
    water_mask = water;
    water_mask(logical(mask)) = 0;
    
    waternum = unique(water_mask);
    waternum(waternum==0) = [];
    
    for w = 1:length(waternum)
        disp(['parcel #' num2str(waternum(w))])
        waterind = find(water==waternum(w));
        temp = zeros(32492,1);
        temp(waterind) = w+1000*(h-1);
        %save(gifti(single(temp)),[dir '/parcel' num2str(w) '.func.gii'])
        %evalc(['!mv ' dir '/parcel' num2str(w) '.func.gii ' dir '/parcel' num2str(w) '.metric']);
        save(gifti(single(temp)),[dir '/parcel.func.gii'])
        evalc(['!mv ' dir '/parcel.func.gii ' dir '/parcel.metric']);
        
        evalc(['!caret_command64 -surface-to-volume ' surfdir '/Conte69.' HEMS{h} '.midthickness.32k_fs_LR.coord.gii ' surfdir '/Conte69.' HEMS{h} '.32k_fs_LR.topo.gii ' dir '/parcel.metric 1 ' dir '/empty_MNI_111.nii']);
        %evalc(['!nifti_4dfp -4 ' dir '/empty_MNI_222.nii ' dir '/empty_MNI_222']);
        %image_temp = read_4dfpimg_HCP([dir '/empty_MNI_222.4dfp.img']);
        data = load_nii([dir '/empty_MNI_111.nii']); image_temp = data.img;
        image_temp(image_temp>0) = waternum(w);%w+1000*(h-1);
        parcels{h}(:,:,:,w) = image_temp;
        
        %community = mean(communities(waterind));
        %dlmwrite(outtextname,[w+1000*(h-1) community],'delimiter',' ','-append')
    end
end


add_parcels = sum(parcels{1},4) + sum(parcels{2},4);
data.img = add_parcels;
save_nii(data,'Parcels_MNI_111.nii')
system('nifti_4dfp -4 Parcels_MNI_111.nii Parcels_MNI_111.4dfp.img')
system('t4img_4dfp none Parcels_MNI_111.4dfp.img Parcels_MNI_222.4dfp.img -n -Oempty_MNI_222.4dfp.ifh')
system('t4img_4dfp none Parcels_MNI_111.4dfp.img Parcels_MNI_333.4dfp.img -n -OMNI152_3mm_brain.4dfp.ifh')

system('t4img_4dfp /data/nil-bluearc/raichle/lin64-tools/MNI152lin_T1_to_711-2B_t4 Parcels_MNI_111.4dfp.img Parcels_711-2b_111.4dfp.img -n -O111')
system('t4img_4dfp /data/nil-bluearc/raichle/lin64-tools/MNI152lin_T1_to_711-2B_t4 Parcels_MNI_111.4dfp.img Parcels_711-2b_222.4dfp.img -n -O222')
system('t4img_4dfp /data/nil-bluearc/raichle/lin64-tools/MNI152lin_T1_to_711-2B_t4 Parcels_MNI_111.4dfp.img Parcels_711-2b_333.4dfp.img -n -O333')

system('nifti_4dfp -n Parcels_MNI_222.4dfp.img Parcels_MNI_222.nii')
system('nifti_4dfp -n Parcels_MNI_333.4dfp.img Parcels_MNI_333.nii')
system('nifti_4dfp -n Parcels_711-2b_111.4dfp.img Parcels_711-2b_111.nii')
system('nifti_4dfp -n Parcels_711-2b_222.4dfp.img Parcels_711-2b_222.nii')
system('nifti_4dfp -n Parcels_711-2b_333.4dfp.img Parcels_711-2b_333.nii')

%write_4dfpimg(add_parcels,[dir '/parcels.4dfp.img'],'bigendian')
%write_4dfpifh_HCP([dir '/parcels.4dfp.img'],1,'bigendian',91,109,91,2,2,2);
