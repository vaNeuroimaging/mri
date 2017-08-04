warning off

subjects = {'396'};
    %'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','181','182'};
%'166'

seedname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/Cluster/pDMN_Prec_roi.mat';

seedimagename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/pDMN_Prec_Clus_funcspace.nii';

%seedimage_nooutlinename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/TNN-aDMN-vmPFC_12mm_4space.nii';
%

brainmaskname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.hdr';


seed = maroi('load_cell',seedname);
brainmask = load_nii(brainmaskname);

FCoutputimgname = 'Prec_FC.nii';
SCoutputimgname = 'Prec_SC.nii.gz';

%CorrelOutputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Func_Struct/vmPFC.nii';

for subject = 1:length(subjects)
    subjid = subjects{subject};
    
    
    data_location = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/SPM8/Rest/'];
    datafile4D =  [data_location 'filt_smoothed.nii.gz'];
    datafiles3D = [data_location 'fsw*.nii'];
    outputfile = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/'];
    
    DTIdir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/'];
    
    try rmdir([outputfile 'temp'],'s');catch;end
    try rmdir([outputfile 'voxeltemp'],'s');catch;end
    
    
%     motionparamfile = dir([data_location 'rp*.txt']);
%     motionparams = textread([data_location motionparamfile.name]);
%     
    datafilenames3D = dir(datafiles3D);
    m = size(datafilenames3D, 1);
    for i=1:m
        P(i, :) = [data_location datafilenames3D(i).name];
    end
    
    
    %VOXELWISE FUNCTIONAL CONNECTIVITY
    %----------------------------------------------------------------------
    
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
    
    
    disp(['Subject ' subjid ': Functional Connectivity'])
    
    %R = jacket_partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams]);
    R = partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams]);
    Fishervals = .5*(log(1+R)-log(1-R));
    Fisherimg(:,:,:) = reshape(Fishervals,size(data.img,1),size(data.img,2),size(data.img,3));
    clear R Fishervals
    
    Fisherimg(find(brainmask.img==0)) = 0;
    
    try mkdir(outputfile); catch; end
    
    try delete([outputfile FCoutputimgname]); catch; end
    
    outputimage = make_nii(Fisherimg,[2 2 2],[40 57 26]);
    save_nii(outputimage,[outputfile FCoutputimgname]);
    
    
    
    
    
    
    clear Fisherimg data reshapeddata motionparams seed_timecourse WM_timecourse CSF_timecourse Fisheroutput string m
    
    
  
    
    %VOXELWISE STRUCTURAL CONNECTIVITY
    %----------------------------------------------------------------------
    
    
    mkdir([outputfile 'temp']);
    delete([outputfile 'temp/masks.txt']);
    fid = fopen([outputfile 'temp/masks.txt'],'at');
    fprintf(fid,'%s\n\r\%s',seedimagename,[outputfile 'temp/temp_target.nii']);
    fclose(fid);
    
    disp(['Subject ' subjid ': Structural Connectivity'])
    
    eval(['!flirt -in ' P(1, :) ' -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -omat ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
    
    eval(['!convert_xfm -omat ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -inverse ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat']);
    
    evalc(['!/data/apps/fsl/4.1.7/bin/probtrackx --mode=seedmask -x ' seedimagename ' -l -c 0.3 -S 4000 --steplength=0.1 -P 5000 --xfm=' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat --forcedir --opd --pd -s ' DTIdir 'bothruns.bedpostX/merged -m ' DTIdir '/bothruns.bedpostX/nodif_brain_mask  --dir=' outputfile 'temp'])

    copyfile([outputfile 'temp/fdt_paths.nii.gz'],[outputfile SCoutputimgname]);
    gunzip([outputfile SCoutputimgname]);
    
    try rmdir([outputfile 'temp'],'s');catch;end
    
    disp(' ');
    
    clear P
       
end


% 
% CorrelOutputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Func_Struct/correl_Prec.nii';
% TOutputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Func_Struct/T_Prec.nii';
% 
% 
% 
% for subject = 1:length(subjects)
%     subjid = subjects{subject};
%     
%     try
%     tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' FCoutputimgname]);
%     catch
%         gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' FCoutputimgname '.gz']);
%         tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' FCoutputimgname]);
%     end
%     
%     fcdata(subject,:,:,:) = tempdata.img;
%     
%     try
%         tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' SCoutputimgname]);
%     catch
%         gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' SCoutputimgname '.gz']);
%         tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' SCoutputimgname]);
%     end
%     
%     scdata(subject,:,:,:) = tempdata.img;
%     
% end
% 
% CorrelOutput = tempdata;
% 
% 
% reshapedfcdata = reshape(fcdata,[length(subjects) size(fcdata,2)*size(fcdata,3)*size(fcdata,4)]);
% reshapedscdata = reshape(scdata,[length(subjects) size(scdata,2)*size(scdata,3)*size(scdata,4)]);
% 
% 
% for vox = 1:size(reshapedfcdata,2)
%     if length(find(reshapedscdata(:,vox)>0)) > (length(subjects)/5);
%         substouse = find(reshapedscdata(:,vox)>0);
%         [rho, p] = corr(reshapedfcdata(substouse,vox),reshapedscdata(substouse,vox));
%         correlation(vox) = p;
%         direction(vox) = sign(rho);
%     else
%         correlation(vox) = 1;
%         direction(vox) = 0;
%     end
% end
% 
% correlationoutput = (1-reshape(correlation,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]))  .* reshape(direction,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]);
% CorrelOutput.img = correlationoutput;
% 
% save_nii(CorrelOutput,CorrelOutputfilename);
% 
% clear direction
% 
% 
% TOutput = tempdata;
% 
% for vox = 1:size(reshapedfcdata,2)
%     if length(find(reshapedscdata(:,vox)>0)) > 2  && length(find(reshapedscdata(:,vox)>0)) < (length(subjects)-2);
%         
%         [H, p, CI, stats] = ttest2(reshapedfcdata(find(reshapedscdata(:,vox)>0),vox), reshapedfcdata(find(reshapedscdata(:,vox)==0),vox));
%         
%         ttest(vox) = p;
%         direction(vox) = sign(stats.tstat);
%     else
%         ttest(vox) = 1;
%         direction(vox) = 0;
%     end
% end
% 
% ttestoutput = (1-reshape(ttest,[size(fcdata,2) size(fcdata,3) size(fcdata,4)])) .* reshape(direction,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]);
% TOutput.img = ttestoutput;
% 
% save_nii(TOutput,TOutputfilename);

