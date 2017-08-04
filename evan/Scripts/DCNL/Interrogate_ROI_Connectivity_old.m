%Interrogate_ROI
%
%Extracts the timecourses within defined ROIs for the subjects.  Regresses
%out the signal from White Matter and from CSF.  Correlates residual
%timecourses against each other.  Converts the Pearson's R values to
%Fisher's Z.  Saves all Z values in a file which is easily exportable into
%Excel.
%
%ROIs and subjects are specified at the top of the script; data location is
%specified in the middle
%
%
%Created by E. Gordon 5/08/08


warning off


%USER INPUT
%--------------------------------------------------------------------------
%Name of output file
outputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Taren et al ROIs/ConnectivityOutput.txt';

%Names of subjects
subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};

%Names of runs to interrogate
runs = {'Rest','Nback'};

%Header of functional data to interrogate
header = 'sw';

%Names of ROIs to interrogate
roi_names = {'DL1','DL2','DL3','DL4','DL5','DL6','DM1','DM2','DM3','DM4','DM5'};

%Correlations between ROIs to run.  Each bracketed number pair is a correlation; each number within the brackets is an ROI (the number reflects the order within roi_names)
comparisons = {[1 7] [1 8] [1 9] [1 10] [1 11] [2 7] [2 8] [2 9] [2 10] [2 11] [3 7] [3 8] [3 9] [3 10] [3 11] [4 7] [4 8] [4 9] [4 10] [4 11] [5 7] [5 8] [5 9] [5 10] [5 11] [6 7] [6 8] [6 9] [6 10] [6 11]};

%Timepoints within functional run to use in correlation.  If all, leave blank.
timepointstouse = [];

%location of ROIs
analysis_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Taren et al ROIs/'];

%END USER INPUT
%--------------------------------------------------------------------------

%Set up output file
delete([outputfilename]);
fid = fopen([ outputfilename],'at');
fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Subject','Run','ROI_Comparison','Correlation');
fclose(fid);
dlmwrite([ outputfilename],' ','-append');


conditions = {'1Back','2Back','3Back'};
%condtimes = {[49 93 137] [5 115 159] [27 71 181]};
condtimes = {[48 92 136] [4 114 158] [26 70 180]};

%%%%%%%%%%%%%%%conddur = 11;
conddur = 12;

alltimepointstouse = [];
for cond = 1:length(conditions)
    timepointstouse{cond} = [condtimes{cond}(1):condtimes{cond}(1)+conddur-1 , condtimes{cond}(2):condtimes{cond}(2)+conddur-1, condtimes{cond}(3):condtimes{cond}(3)+conddur-1];
    alltimepointstouse = [alltimepointstouse timepointstouse{cond}];
end


load('taskregressor.mat');
taskregressor = taskregressor(alltimepointstouse,:);



%Subject loop
for subject = 1:length(subjects)

    clear timecourse residual_timecourse residuals

    clear datamatrix;
    
    %Run loop
    for run = 1:length(runs)

        subjsid = subjects{subject};
        
        %USER INPUT
        %--------------------------------------------------------------------------
        %location of functional data
        data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/' runs{run} '/'];

        %Location of CSF and WM ROIs
        CSFROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_CSF_roi.mat'];
        WMROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_WM_roi.mat'];
        
        %END USER INPUT
        %--------------------------------------------------------------------------

        motionparamfile = dir([data_dir 'rp*.txt']);
        motionparams = textread([data_dir motionparamfile.name]);
        
        clear P Y
        
        %Find functional images
        imgfiles = dir([data_dir header '*.hdr']);
        m = size(imgfiles, 1);
        for j=1:m
            P(j, :) = [data_dir imgfiles(j).name];
        end

        if strcmp(runs{run},'Nback')
            P = P(alltimepointstouse,:);
            motionparams = motionparams(alltimepointstouse,:);
        end

        %Load and interrogate CSF ROI
        disp([subjsid ' : ' runs{run} ' : CSF'])
        LV_rois = maroi('load_cell',CSFROI);
        [Y a b c] = getdata(LV_rois{1}, P,'l');
        LV_timecourse{run} = mean(Y,2);
        clear Y

        %Load and interrogate WM ROI
        disp([subjsid ' : ' runs{run} ' : WM'])
        WM_rois = maroi('load_cell',WMROI);
        [Y a b c] = getdata(WM_rois{1}, P,'l');
        WM_timecourse{run} = mean(Y,2);
        clear Y

        %ROI loop
        for roinum = 1:length(roi_names)
            disp([subjsid ' : ' runs{run} ' : ' roi_names{roinum}])
            
            %Load and interrogate ROI of interest
            roi = [analysis_dir, roi_names{roinum} , '_roi.mat'];
            rois = maroi('load_cell', roi);
            [Y a b c] = getdata(rois{1}, P,'l');
            timecourse{run,roinum} = mean(Y,2);

            %Regress CSF and WM signal out of every voxel
            for voxelnum = 1:size(Y,2)
                if strcmp(runs{run},'Nback')
                    [b, bint, r] = regress(Y(:,voxelnum),[LV_timecourse{run} WM_timecourse{run} motionparams taskregressor ones(length(LV_timecourse{run}),1)]);
                else
                    [b, bint, r] = regress(Y(:,voxelnum),[LV_timecourse{run} WM_timecourse{run} motionparams ones(length(LV_timecourse{run}),1)]);
                end
                residuals(:,voxelnum) = r;
            end
            
            residual_timepoint_means = mean(residuals,2);
            clear residuals filtered_residuals

            %Average across all voxels
            residual_timecourse{run,roinum} = residual_timepoint_means;

        end
    end
    
    
    for comparison = 1:length(comparisons)

        for run = 1:length(runs)

            %Correlate residual timecourses, convert R to Z, and write to output file
            correlvals = corrcoef(residual_timecourse{run,comparisons{comparison}(1)},residual_timecourse{run,comparisons{comparison}(2)});
            Fisher = .5*(log(1+correlvals(2,1))-log(1-correlvals(2,1)));
            texttowrite = [subjsid,'   ' runs{run} '   ',[roi_names{comparisons{comparison}(1)} '_vs_' roi_names{comparisons{comparison}(2)}],'   ',num2str(Fisher)];
            dlmwrite([ outputfilename],texttowrite,'-append','delimiter','');

        end

    end
end
