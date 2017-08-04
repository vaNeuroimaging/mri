subjects = {'vc34096' 'vc34116' 'vc34125' 'vc34126' 'vc34128' 'vc34140' 'vc34141' 'vc34198' 'vc34199' 'vc34200' 'vc34201' 'vc34220' 'vc34250' 'vc34251' 'vc34252' 'vc34253' 'vc34268' 'vc34306' 'vc34307' 'vc34308' 'vc34330' 'vc34331'  'vc34401' 'vc34402' 'vc34403' 'vc34404' 'vc34405' 'vc34408'};
%
for subject = 1:length(subjects)
    try %evalc(['!cp /data/cn3/joe/ResourceDataLimited/' subjects{subject} '/atlas/' subjects{subject} '_mpr_n1_111_t88.4dfp.* ./']);
        
        %evalc(['!nifti_4dfp -n ' subjects{subject} '_mpr_n1_111_t88.4dfp.img ' subjects{subject} '_mpr_n1_111_t88.nii']);
        %copyfile([subjects{subject} '_mpr_n1_111_t88.nii'],['/data/cn4/segmentation/freesurfer5_supercomputer/' subjects{subject} '/' subjects{subject} '_mpr_n1_111_t88.nii']);
        
        disp(subjects{subject})
        
        evalc(['!fslchfiletype NIFTI_GZ /data/cn4/segmentation/freesurfer5_supercomputer/' subjects{subject} '/' subjects{subject} '_mpr_n1_111_t88.nii']);
        
        if exist(['/data/cn4/segmentation/freesurfer5_supercomputer/' subjects{subject} '/' subjects{subject} '_mpr_n1_111_t88.nii.gz'])
            delete(['/data/cn4/segmentation/freesurfer5_supercomputer/' subjects{subject} '/' subjects{subject} '_mpr_n1_111_t88.nii']);
        end
        
    catch
        disp([subjects{subject} ' failure']);
    end
end