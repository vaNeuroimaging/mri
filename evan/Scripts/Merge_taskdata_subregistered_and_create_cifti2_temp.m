outputdir = '/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/Katie/';


lists = {'/data/cn4/katie/DEV_PHONOLOGY_August2010/REPREPROCESS/GLM_lists/glmv1R_32adults.glm_list';};

directories = {'/data/cn4/evan/Task_parcellation/Ourdata_unsmoothed/Katie/SingleSubjects/'};%,...

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';

hems = {'L' 'R'};


Mergedata{1} = [];
Mergedata{2} = [];
for listnum = 1:length(lists)
    fslmergestr = 'fslmerge -t Manytasks.nii.gz Manytasks.nii.gz';
    list = lists{listnum};
    disp(list)
    directory = directories{listnum};
    
    [front entries] = textread(list,'%s%s','delimiter',':');
    
    %Figure out the number of subjects, the glm files, and the T4 files from
    %the GLM list file
    numsubjects = str2num(entries{2});
    glmfiles = entries(3:2+numsubjects);
    T4files = entries(4+numsubjects:end);
    
    %This part accounts for different GLM list styles
    for i = 1:length(glmfiles)
        if strcmp(glmfiles{i},'')
            glmfiles{i} = front{2+i};
        end
        if strcmp(T4files{i},'')
            T4files{i} = front{3+numsubjects+i};
        end
    end
    
    for subject = 1:length(glmfiles)
        disp(glmfiles{subject})
        cd(outputdir)
        subnameindex = strfind(glmfiles{subject},'/vc');
        
        slashindex = min([strfind(glmfiles{subject}(subnameindex(1)+1:end),'/') strfind(glmfiles{subject}(subnameindex(1)+1:end),'_') strfind(glmfiles{subject}(subnameindex(1)+1:end),'.')]);
        
        subname{subject} = glmfiles{subject}(subnameindex+1:subnameindex+slashindex-1);
        
        surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR';
        HEMS = {'L';'R'};
        smooth = '2.55';
        
        conditionfiles = dir([directory subname{subject} '*.4dfp.img']);
        
        submask = ['/data/cn4/evan/Task_parcellation/Manytasks_subregister/goodvoxels/' subname{subject} '_goodvoxels.nii.gz'];
        
        for condition = 1:length(conditionfiles)
            
            subfunc = [directory '/' conditionfiles(condition).name(1:end-9)];
            system(['niftigz_4dfp -n ' subfunc ' ' subfunc])
            
            for hem = 1:2
                
%                 if hem==1
%                     
%                     if listnum == 1 && subject == 1 && condition == 1
%                         copyfile([directory conditionfiles(condition).name(1:end-9) '.nii.gz'],'Manytasks.nii.gz')
%                     else
%                         fslmergestr = [fslmergestr ' ' directory conditionfiles(condition).name(1:end-9) '.nii.gz'];
%                     end
%                 end
                
                midsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.midthickness.native.surf.gii'];
                midsurf_LR32k = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
                whitesurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.white.native.surf.gii'];
                pialsurf = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.pial.native.surf.gii'];
                nativedefsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/Native/' subname{subject} '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
                outsphere = [surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/' subname{subject} '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        
                surfname = [subname{subject} '_' HEMS{hem} '_condition' num2str(condition)];
                disp('Map volume to surface')
                system([workbenchdir '/wb_command -volume-to-surface-mapping ' subfunc '.nii.gz ' midsurf ' ' outputdir '/' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -volume-roi ' submask]);
                disp('Dilate surface timecourse')
                system([workbenchdir '/wb_command -metric-dilate ' outputdir '/' surfname '.func.gii ' midsurf ' 10 ' outputdir '/' surfname '_dil10.func.gii'])
                
                disp('Deform timecourse to 32k fs_LR')
                cd([surfdir '/' subname{subject} '/7112b_fs_LR/fsaverage_LR32k/'])
                %system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' outputdir '/' surfname '_dil10.func.gii ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii']);
                system([workbenchdir '/wb_command -metric-resample ' outputdir '/' surfname '_dil10.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
                
                disp('Smooth surface timecourse')
                system([workbenchdir '/wb_command -metric-smoothing ' midsurf_LR32k ' ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii ' smooth ' ' outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii'])
                
                %thisdata = gifti([outputdir '/' surfname '_dil10_32k_fsLR_smooth' smooth '.func.gii']);
                %Mergedata{hem} = [Mergedata{hem} thisdata.cdata];
                
                system(['rm ' outputdir '/' surfname '.func.gii'])
                system(['rm ' outputdir '/' surfname '_dil10.func.gii'])
                system(['rm ' outputdir '/' surfname '_dil10_32k_fsLR.func.gii'])
            end
        end
        
        
        
    end
    clear subname
    system(fslmergestr)
end

%save(gifti(single(Mergedata{1})),'Manytasks_L.func.gii','ExternalFileBinary')
%save(gifti(single(Mergedata{2})),'Manytasks_R.func.gii','ExternalFileBinary')
%
%gunzip Manytasks.nii.gz
%system('wb_command -cifti-create-dense-timeseries Manytasks.dtseries.nii -volume Manytasks.nii /data/cn4/laumannt/subcortical_mask/mode_subcortical_label_LR_333.nii -left-metric Manytasks_L.func.gii -roi-left /data/cn4/evan/RestingState/FC_Mapping_120_nocrap/L.atlasroi_erode3_nocrap.32k_fs_LR.shape.gii -right-metric Manytasks_R.func.gii -roi-right /data/cn4/evan/RestingState/FC_Mapping_120_nocrap/R.atlasroi_erode3_nocrap.32k_fs_LR.shape.gii')
