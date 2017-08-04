warning off

subjects = {'400','401','402','406','407','410','412','415','416','417','420'};
%{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295'};
%'297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
    %
    %
%


seedname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/TNN-aDMN-vmPFC_12mm_roi.mat';

seedimagename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/TNN-aDMN-vmPFC_12mm_4space_outline.nii';

seedimage_nooutlinename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/TNN-aDMN-vmPFC_12mm_4space.nii';
%

brainmaskname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.hdr';


seed = maroi('load_cell',seedname);
brainmask = load_nii(brainmaskname);

FCoutputimgname = 'vmPFC_FC.nii';
SCoutputimgname = 'vmPFC_SC.nii';

%CorrelOutputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Func_Struct/vmPFC.nii';

for subject = 1:length(subjects)
    subjid = subjects{subject};
    
    
    data_location = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/SPM8/FirstRest/'];
    datafile4D =  [data_location 'filt_smoothed.nii.gz'];
    datafiles3D = [data_location 'fsw*'];
    outputfile = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/'];
    
    DTIdir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/'];
    
    try rmdir([outputfile 'temp'],'s');catch;end
    try rmdir([outputfile 'voxeltemp'],'s');catch;end
    
    disp(['Subject ' subjid ': Functional Connectivity'])
    
    motionparamfile = dir([data_location 'rp*.txt']);
    motionparams = textread([data_location motionparamfile.name]);
    
    datafilenames3D = dir([datafiles3D '.hdr']);
    if isempty(datafilenames3D)
        datafilenames3D = dir([datafiles3D '.nii']);
    end
    m = size(datafilenames3D, 1);
    for i=1:m
        P(i, :) = [data_location datafilenames3D(i).name];
    end
    
    [Y a b c] = getdata(seed{1}, P,'l');
    seed_timecourse = mean(Y,2);
    
    CSF_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjid '_CSF_roi.mat']);
    [Y a b c] = getdata(CSF_rois{1}, P,'l');
    CSF_timecourse = mean(Y,2);
    
    WM_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjid '_WM_roi.mat']);
    [Y a b c] = getdata(WM_rois{1}, P,'l');
    WM_timecourse = mean(Y,2);
    
    try gunzip(datafile4D); catch; end;
    
    data = load_nii(datafile4D(1:end-3));
    
    reshapeddata = double(reshape(data.img,size(data.img,1)*size(data.img,2)*size(data.img,3),size(data.img,4))');
    
    
    %R = jacket_partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams]);
    R = partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams]);
    Fishervals = .5*(log(1+R)-log(1-R));
    Fisherimg(:,:,:) = reshape(Fishervals,size(data.img,1),size(data.img,2),size(data.img,3));
    clear R Fishervals
    
    Fisherimg(find(brainmask.img==0)) = 0;
    
    try mkdir(outputfile); catch; end
    
    try delete([outputfile FCoutputimgname]); catch; end
    
    outputimage = make_nii(Fisherimg,data.hdr.dime.pixdim(2:4), data.hdr.hist.originator(1:3));
    save_nii(outputimage,[outputfile FCoutputimgname]);
    
    clear Fisherimg data reshapeddata motionparams seed_timecourse WM_timecourse CSF_timecourse Fisheroutput string m
    
    
  
     mkdir([outputfile 'temp']);
    delete([outputfile 'temp/masks.txt']);
    fid = fopen([outputfile 'temp/masks.txt'],'at');
    fprintf(fid,'%s\n\r\%s',seedimagename,[outputfile 'temp/temp_target.nii']);
    fclose(fid);
    mkdir([outputfile 'voxeltemp']);
    
    disp(['Subject ' subjid ': First-pass Structural Connectivity'])
    
%     evalc(['!flirt -in /fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/mean_FA_skeleton_mask_inv_thr4.nii.gz -ref /fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/origdata/' subjid '_dti_FA.nii.gz -out ' DTIdir 'bothruns/transformed_inverse_skeleton -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
%     evalc(['!fslmaths ' DTIdir 'bothruns/transformed_inverse_skeleton -thr .5 -bin ' DTIdir 'bothruns/transformed_inverse_skeleton']);
    
    eval(['!flirt -in ' DTIdir 'bothruns/nodif_brain.nii.gz -ref /fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/MNI152_T1_4mm_brain_converted.nii -omat ' DTIdir 'bothruns.bedpostX/xfms/diff24mm.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
%     evalc(['!flirt -in /fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/mean_FA_skeleton_mask_inv_thr4.nii.gz -ref /fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/MNI152_T1_4mm_brain_converted.nii -out /fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/mean_FA_skeleton_mask_inv_thr4_4mm.nii.gz -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
%     evalc(['!fslmaths /fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/mean_FA_skeleton_mask_inv_thr4_4mm.nii.gz -thr .5 -bin /fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/mean_FA_skeleton_mask_inv_thr4_4mm.nii.gz']);
    
    eval(['!convert_xfm -omat ' DTIdir 'bothruns.bedpostX/xfms/4mm2diff.mat -inverse ' DTIdir 'bothruns.bedpostX/xfms/diff24mm.mat']);
    
    eval(['!flirt -in ' DTIdir '/bothruns/dti_FA.nii.gz -ref /fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/MNI152_T1_4mm_brain_converted.nii  -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/diff24mm.mat -out ' outputfile 'temp/FAmap']);
   %eval(['!flirt -in ' DTIdir '/bothruns/dti_FA.nii.gz -ref ' P(1, :) ' -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -out ' outputfile 'temp/FAmap']);
    
    evalc(['!/data/apps/fsl/4.1.7/bin/probtrackx --mode=seedmask -x ' seedimagename ' -l -c 0.2 -S 1000 --steplength=0.5 -P 5000 --xfm=' DTIdir 'bothruns.bedpostX/xfms/4mm2diff.mat --forcedir --opd -s ' DTIdir 'bothruns.bedpostX/merged -m ' DTIdir '/bothruns.bedpostX/nodif_brain_mask  --dir=' outputfile 'temp']);
    
    eval(['!fslmaths ' outputfile 'temp/fdt_paths.nii.gz -mul ' seedimage_nooutlinename ' ' outputfile 'temp/fdt_paths_within_mask.nii.gz']);
    eval(['!fslmaths ' outputfile 'temp/fdt_paths.nii.gz -sub ' outputfile 'temp/fdt_paths_within_mask.nii.gz ' outputfile 'temp/fdt_paths_outside_mask.nii.gz']);
    
    eval(['!fslchfiletype NIFTI ' outputfile 'temp/fdt_paths_outside_mask.nii.gz']);
    
    Paths = load_nii([outputfile 'temp/fdt_paths_outside_mask.nii']);
    
    pathpresentindices = find(Paths.img);
    Onevoxel = Paths;
    StructConnectData = double(zeros(size(Onevoxel.img)));
    
    for voxel = 1:length(pathpresentindices);
        
        string{voxel} = ['Subject ' subjid ': Voxelwise Structural Connectivity: voxel number ' num2str(voxel) ' of ' num2str(length(pathpresentindices))];
        if voxel==1; fprintf('%s',string{voxel}); else; fprintf([repmat('\b',1,length(string{voxel-1})) '%s'],string{voxel}); end
        
        Onevoxel.img = zeros(size(Onevoxel.img));
        Onevoxel.img(pathpresentindices(voxel)) = 1;
        save_nii(Onevoxel,[outputfile 'temp/temp_target.nii']);
        evalc(['!/data/apps/fsl/4.1.7/bin/probtrackx --network --mode=seedmask -x ' outputfile 'temp/temp_target.nii --waypoints=' seedimagename ' --stop=' seedimagename ' -l -c 0.2 -S 500 --steplength=1 -P 250 --xfm=' DTIdir 'bothruns.bedpostX/xfms/4mm2diff.mat --forcedir --opd -s ' DTIdir 'bothruns.bedpostX/merged -m ' DTIdir 'bothruns.bedpostX/nodif_brain_mask  --dir=' outputfile 'voxeltemp']);
        evalc(['!fslmaths ' outputfile 'voxeltemp/fdt_paths.nii.gz -bin -mul ' outputfile 'temp/FAmap.nii.gz ' outputfile 'voxeltemp/FA_within_tract']);
        meanFA = evalc(['!fslstats ' outputfile 'voxeltemp/FA_within_tract.nii.gz -M']);
        
        if meanFA == 0
            evalc(['!/data/apps/fsl/4.1.7/bin/probtrackx --network --mode=seedmask -x ' outputfile 'temp/temp_target.nii --waypoints=' seedimagename ' --stop=' seedimagename ' -l -c 0.2 -S 1000 --steplength=.5 -P 1000 --xfm=' DTIdir 'bothruns.bedpostX/xfms/4mm2diff.mat --forcedir --opd -s ' DTIdir 'bothruns.bedpostX/merged -m ' DTIdir 'bothruns.bedpostX/nodif_brain_mask  --dir=' outputfile 'voxeltemp']);
            evalc(['!fslmaths ' outputfile 'voxeltemp/fdt_paths.nii.gz -bin -mul ' outputfile 'temp/FAmap.nii.gz ' outputfile 'voxeltemp/FA_within_tract']);
            meanFA = evalc(['!fslstats ' outputfile 'voxeltemp/FA_within_tract.nii.gz -M']);
        end
        
        StructConnectData(pathpresentindices(voxel)) = str2num(meanFA);
    
    end
        
    StructConnectOutput = make_nii(StructConnectData,Paths.hdr.dime.pixdim(2:4),Paths.hdr.hist.originator(1:3),16);
    save_nii(StructConnectOutput,[outputfile '4mmspace_' SCoutputimgname]);
    
    
    eval(['!flirt -ref /fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/MNI152_T1_4mm_brain_converted.nii -in ' P(1, :) ' -omat ' outputfile 'temp/func24mm -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
    eval(['!flirt -in ' outputfile FCoutputimgname ' -ref /fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/MNI152_T1_4mm_brain_converted.nii -applyxfm -init ' outputfile 'temp/func24mm -out ' outputfile '4mmspace_' FCoutputimgname]);
    eval(['!fslchfiletype NIFTI ' outputfile '4mmspace_' SCoutputimgname]);
    
    try rmdir([outputfile 'temp'],'s');catch;end
    try rmdir([outputfile 'voxeltemp'],'s');catch;end
    
    disp(' ');
    
    clear P
       
end
