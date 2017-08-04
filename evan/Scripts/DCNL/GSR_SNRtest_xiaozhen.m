subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};

mask_name='/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.hdr';

mask=load_nii(mask_name);

mask=mask.img;
mask=single(mask);
mask_ind=find(mask==1); % the index within the mask


non_zero_l=length(mask_ind);
gni =zeros(length(subjects),length(mask_ind),2);
gni_sum=zeros(length(subjects),1);


for fi=1:length(subjects)
    
    subname= subjects{fi};
    disp(subname)
    
    data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subname, '/SPM8/FirstRest/'];
    
    
    smdata=[data_dir 'filt_smoothed.nii'];
    datav=load_nii(smdata);
    datav=datav.img;
    
    imgfiles = dir([data_dir 'fsw*.nii']);
    %Put image names into a char array
    m = size(imgfiles, 1);
    for j=1:m
        P(j, :) = [data_dir imgfiles(j).name];
    end
    
    
    CSFROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subname '_CSF_roi.mat'];
    WMROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subname '_WM_roi.mat'];
    
    donegetdata = 0;
    
    %Be robust to getdata crashing (which it does a lot recently)
    while donegetdata == 0
        try
            
            %Load and interrogate CSF ROI
            CSF_rois = maroi('load_cell',CSFROI);
            [Y a b c] = getdata(CSF_rois{1}, P,'l');
            CSF_timecourse = mean(Y,2);
            clear Y
            
            %Load and interrogate WM ROI
            WM_rois = maroi('load_cell',WMROI);
            [Y a b c] = getdata(WM_rois{1}, P,'l');
            WM_timecourse = mean(Y,2);
            clear Y
            
            donegetdata = 1;
            
        catch
            
            disp(['Getdata failed: ' datestr(now)])
            
        end
    end
    
    %Find and load motion param file
    motionparamfile = dir([data_dir 'rp*.txt']);
    motionparams = textread([data_dir motionparamfile.name]);
    
    %Add the calculated WM/CSF timecourses as regressors
    regressors = [CSF_timecourse WM_timecourse motionparams];
    
    
    time_point = size(datav,4);
    
    data=zeros(non_zero_l,time_point,'single');
    
    for ti=1:time_point
        volume_i=datav(:,:,:,ti);
        data(:,ti)=volume_i(mask_ind);
    end
    glob_sig =mean(data,1);
    data=data';
    glob_sig=glob_sig';
    
    for voxi=1:length(mask_ind)
        
        [r,p] = partialcorr(data(:,voxi),glob_sig, regressors);
        gni(fi,voxi,1) = r;
        gni(fi,voxi,2) =p;
    end
    gni_sum(fi)=length(find(gni(fi,:,1)<0 & gni(fi,:,2)<0.05))/length(mask_ind);
    
    clear data datav glob_sig volume_i P regressors
    
end




figure;

plot(1:length(subjects),gni_sum,'+');
%save('asd_task_sub_gni.mat', 'gni_sum', 'gni');


