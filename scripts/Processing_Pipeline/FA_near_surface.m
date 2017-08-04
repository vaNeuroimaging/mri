function FA_near_surface(subject,distances)


Lpial = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.pial.32k_fs_LR.surf.gii']);
Rpial = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.pial.32k_fs_LR.surf.gii']);
Lwhite = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.white.32k_fs_LR.surf.gii']);
Rwhite = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.white.32k_fs_LR.surf.gii']);

for d = 1:length(distances)
    
    distance = distances(d);
    
    Linteriorsurf = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.interior' num2str(distance) '.32k_fs_LR.surf.gii'];
    Rinteriorsurf = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.interior' num2str(distance) '.32k_fs_LR.surf.gii'];
    
    if ~exist(Linteriorsurf,'file') || ~exist(Rinteriorsurf,'file')
        
        vecsL = Lwhite.vertices - Lpial.vertices;
        mags = sqrt(sum((vecsL.^2),2));
        normvecsL = vecsL .* distance ./ repmat(mags,1,3);
        interiorL = Lwhite.vertices + normvecsL;
        
        vecsR = Rwhite.vertices - Rpial.vertices;
        mags = sqrt(sum((vecsR.^2),2));
        normvecsR = vecsR .* distance ./ repmat(mags,1,3);
        interiorR = Rwhite.vertices + normvecsR;
        
        outL = Lwhite;
        outL.vertices = interiorL;
        save(outL,Linteriorsurf);
        
        outR = Rwhite;
        outR.vertices = interiorR;
        save(outR,Rinteriorsurf);
        
    end
    
    medial_mask_L = '/home/data/scripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii';
    medial_mask_R = '/home/data/scripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';
    
        system(['wb_command -volume-to-surface-mapping /home/data/subjects/' subject '/DTI/DTI_avg_ec_FA_MNI.nii.gz /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.midthickness.32k_fs_LR.surf.gii /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_L.func.gii -ribbon-constrained ' Linteriorsurf ' /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.white.32k_fs_LR.surf.gii']);
        system(['wb_command -volume-to-surface-mapping /home/data/subjects/' subject '/DTI/DTI_avg_ec_FA_MNI.nii.gz /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.midthickness.32k_fs_LR.surf.gii /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_R.func.gii -ribbon-constrained ' Rinteriorsurf ' /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.white.32k_fs_LR.surf.gii']);
        system(['wb_command -cifti-create-dense-timeseries /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_LR.dtseries.nii -left-metric /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_L.func.gii -roi-left ' medial_mask_L ' -right-metric /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_R.func.gii -roi-right ' medial_mask_R]);
        if d==1
            cifti_out = ft_read_cifti_mod(['/home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_LR.dtseries.nii']);
            cifti_out.dimord = 'scalar_pos';
        end
        thisdistdata = ft_read_cifti_mod(['/home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_LR.dtseries.nii']);
        cifti_out.data(:,d) = thisdistdata.data;
        cifti_out.mapname{d} = ['FA within ' num2str(distance) 'mm of white/gray boundary'];
    
end

delete(['/home/data/subjects/' subject '/DTI/' subject '.FA_within_*'])
ft_write_cifti_mod(['/home/data/subjects/' subject '/DTI/FA_near_cortex_LR.dscalar.nii'],cifti_out)
