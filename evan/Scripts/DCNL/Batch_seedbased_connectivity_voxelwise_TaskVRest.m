warning off

subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
%
seedname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_79.gica/ROIs/6mm/TPN-Sal-dACC_roi.mat';

windowsize = 30;

brainmaskname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.hdr';


seed = maroi('load_cell',seedname);
brainmask = load_nii(brainmaskname);


runs = {'Rest','Nback'};
conditions = {'1Back','2Back','3Back'};
condtimes = {[48 92 136] [4 114 158] [26 70 180]};
conddur = 12;

alltimepointstouse = [];
for cond = 1:length(conditions)
    timepointstouse{cond} = [condtimes{cond}(1):condtimes{cond}(1)+conddur-1 , condtimes{cond}(2):condtimes{cond}(2)+conddur-1, condtimes{cond}(3):condtimes{cond}(3)+conddur-1];
    alltimepointstouse = [alltimepointstouse timepointstouse{cond}];
end

load('taskregressor.mat');
taskregressor = taskregressor(alltimepointstouse,:);



for subject = 1:length(subjects)
    subjid = subjects{subject};
    %disp(subjid);
    
    for run = 1:length(runs)
        
        data_location = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/SPM8/' runs{run} '/'];
        
        datafile4D =  [data_location 'smoothed.nii.gz'];
        datafiles3D = [data_location 'sw*.hdr'];
        outputfile = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/VoxelwiseConnectivity/RestVTask/'];
        outputimgname = ['dACC_' runs{run} '.nii'];
        
        motionparamfile = dir([data_location 'rp*.txt']);
        motionparams = textread([data_location motionparamfile.name]);
        
        datafilenames3D = dir(datafiles3D);
        m = size(datafilenames3D, 1);
        for i=1:m
            P(i, :) = [data_location datafilenames3D(i).name];
        end
        
        if strcmp(runs{run},'Nback')
            P = P(alltimepointstouse,:);
            motionparams = motionparams(alltimepointstouse,:);
        else
            P = P(1:108,:);
            motionparams = motionparams(1:108,:);
        end

        getdataworked = 0;
        while getdataworked == 0
            try
                
                [Y a b c] = getdata(seed{1}, P,'l');
                seed_timecourse = mean(Y,2);
                
                CSF_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjid '_CSF_roi.mat']);
                [Y a b c] = getdata(CSF_rois{1}, P,'l');
                CSF_timecourse = mean(Y,2);
                
                WM_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjid '_WM_roi.mat']);
                [Y a b c] = getdata(WM_rois{1}, P,'l');
                WM_timecourse = mean(Y,2);
                
                getdataworked = 1;
                
            catch
               disp('Getdata failure; retrying.') 
            end
        end
        
        try gunzip(datafile4D); catch; end;
        
        data = load_nii(datafile4D(1:end-3));
        datamatrix = data.img;
        
        if strcmp(runs{run},'Nback')
            datamatrix = datamatrix(:,:,:,alltimepointstouse);
        else
            datamatrix = datamatrix(:,:,:,1:108);
        end
        
        
        reshapeddata = double(reshape(datamatrix,size(datamatrix,1)*size(datamatrix,2)*size(datamatrix,3),size(datamatrix,4))');
        
        disp(['Subject ' subjid ', ' runs{run} ': ' outputimgname ])
        
        if strcmp(runs{run},'Nback')
            R = partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams taskregressor]);
        else
            R = partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams]);
        end
        Fishervals = .5*(log(1+R)-log(1-R));
        Fisherimg = reshape(Fishervals,size(datamatrix,1),size(datamatrix,2),size(datamatrix,3));
        clear R Fishervals
            
        
        
        Fisherimg(find(brainmask.img==0)) = 0;
        
        try mkdir(outputfile); catch; end
        
        try delete([outputfile outputimgname]); catch; end
        
        outputimage = make_nii(Fisherimg,[2 2 2],data.hdr.hist.originator(1:3));
        save_nii(outputimage,[outputfile outputimgname]);
        
        clear output data datamatrix reshapeddata motionparams seed_timecourse WM_timecourse CSF_timecourse Fisheroutput string P
        
        
    end
end


