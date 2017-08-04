%Interrogate_ROI
%
%Extracts the timecourses within defined ROIs for the subjects.  Saves those timecourses in in one big file.
%
%ROIs and subjects are specified at the top of the script
%
%
%Created by E. Gordon 5/08/08


warning off

directory = pwd;

outputfilename = 'ConnectivityOutput.txt';

subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274'};


genematrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Data/Genes.txt');

for sub = 1:length(subjects);
    COMTgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),3);
    DATgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),2);
    DRD4genotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),4);
    BDNFgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),5);
    SRTTgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),6);
end

behavior{1} = [484.67	678.21	682.25	568.5	482.8	411.5	853.13	582.5	567.87	555.88	470.94	672.06	430.67	525.13	404.27	488.94	606.64	393.81	792.43	587.47	791	490.88	473.38	559.67	452.23	465.69	634.07	393.93	638	530.94	551.07	578.14	423	921.69	632.94	1040.31	1081.75	460.69];
behavior{2} = [0.93	0.87	1	1	0.93	1	1	1	0.93	1	1	1	0.9	0.92	0.89	1	0.87	1	0.85	0.9	1	1	1	0.93	0.75	1	0.94	0.82	0.93	1	0.92	0.82	0.98	1	0.98	1	1	1];
behavior{3} = [0	0	1	1	0	1	1	1	0	1	1	1	0	0	0	1	0	1	0	0	1	1	1	0	0	1	0	0	0	1	0	0	0	1	0	1	1	1];
behavior{4} = [0.76	0.29	0.41	0.21	0.35	0.35	0.46	0.22	0.26	0.34	0.37	0.55	0.35	0.38	0.16	0.3	0.59	0.43	0.36	0.38	0.56	0.48	0.21	0.2	0.27	0.28	0.57	0.3	0.59	0.41	0.47	0.58	0.3	0.38	0.46	0.54	0.47	0.32];

behavname = {'RT','Acc','Acc_100','ICV'};

runs = {'Nback'};
conditions = {'1Back','2Back','3Back'};
condtimes = {[49 93 137] [5 115 159] [27 71 181]};
conddur = 11;

alltimepointstouse = [];
for cond = 1:length(conditions)
    timepointstouse{cond} = [condtimes{cond}(1):condtimes{cond}(1)+conddur-1 , condtimes{cond}(2):condtimes{cond}(2)+conddur-1, condtimes{cond}(3):condtimes{cond}(3)+conddur-1];
    alltimepointstouse = [alltimepointstouse timepointstouse{cond}];
end


header = 'sw';

baseline_roi_name = 'LV';

analysis_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/ICA_ROIs/'];

roi_names = {'rightFPC_right_dlPFC','rightFPC_right_ilPFC','rightFPC_right_inf_Parietal','leftFPC_left_dlPFC','leftFPC_left_ilPFC','leftFPC_left_inf_parietal','setMaintenance_ACC','setMaintenance_left_insula','setMaintenance_right_insula','dorsalAttention_left_precuneus','dorsalAttention_right_precuneus','anteriorDMN_precuneus','anteriorDMN_vmPFC','posteriorDMN_left_angular','posteriorDMN_precuneus','posteriorDMN_right_angular'};
%'language_left_inf_frontal','language_left_mid_temporal','language_right_inf_frontal','language_right_sup_temporal',,'memory_left_hippocampus','memory_right_hippocampus'
load('Fishermatrix_byload.mat')

for condition = 1:length(conditions)
    
    
%     delete([analysis_dir outputfilename]);
%     fid = fopen([analysis_dir outputfilename],'at');
%     fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Subject','Condition','ROI_Comparison','Fisher');
%     fclose(fid);
%     dlmwrite([analysis_dir outputfilename],' ','-append');
%     for subject = 1:length(subjects)
        
%         clear timecourse residual_timecourse residuals
%         
%         clear datamatrix;
%         
%         subjsid = subjects{subject};
%         
%         data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/' runs{1} '/'];
%         
%         clear roi_files;
%         clear des_path;
%         clear rois;
%         clear des;
%         clear mY;
%         clear ttfile;
%         clear imgfiles;
%         clear n;
%         clear P;
%         clear raw_d;
%         clear fmri_raw_data;
%         clear numVox;
%         
%         imgfiles = dir([data_dir header '*.hdr']);
%         m = size(imgfiles, 1);
%         for j=1:m
%             P(j, :) = [data_dir imgfiles(j).name];
%         end
%         P = P(timepointstouse{condition},:);
%         
%         
%         disp([subjsid ' : ' conditions{condition} ' : LV'])
%         LV_rois = maroi('load_cell',[analysis_dir '../Ind_ROIs/' subjsid '_LV_roi.mat']);
%         [Y a b c] = getdata(LV_rois{1}, P,'l');
%         LV_timecourse{condition} = mean(Y,2);
%         
%         clear Y
%         
%         disp([subjsid ' : ' conditions{condition} ' : WM'])
%         WM_rois = maroi('load_cell',[analysis_dir '../Ind_ROIs/' subjsid '_WM_roi.mat']);
%         [Y a b c] = getdata(WM_rois{1}, P,'l');
%         WM_timecourse{condition} = mean(Y,2);
%         
%         clear Y
%         
%         for roinum = 1:length(roi_names)
%             disp([subjsid ' : ' conditions{condition} ' : ' roi_names{roinum}])
%             roi = [analysis_dir, roi_names{roinum} , '_roi.mat'];
%             %roi = [analysis_dir, roi_names{roinum}];
%             rois = maroi('load_cell', roi);
%             
%             [Y a b c] = getdata(rois{1}, P,'l');
%             
%             timecourse{condition,roinum} = mean(Y,2);
%             
%             for voxelnum = 1:size(Y,2)
%                 [b, bint, r] = regress(Y(:,voxelnum),[LV_timecourse{condition} WM_timecourse{condition} ones(length(LV_timecourse{condition}),1)]);
%                 %[b, bint, r] = regress(Y(:,voxelnum),[WB_timecourse{condition} ones(length(WB_timecourse{condition}),1)]);
%                 residuals(:,voxelnum) = r;
%                 %filtered_residuals(:,voxelnum) = filter(filterparamb,filterparama,r);
%             end
%             
%             residual_timepoint_means = mean(residuals,2);
%             clear residuals filtered_residuals
%             
%             Timecourses{condition}(:,roinum) = residual_timepoint_means;
%             
%         end
%         
%         
%         Correlvals = corrcoef(Timecourses{condition}(:,:));
%         
%         Fishermatrix(:,:,subject,condition) = .5*(log(1+Correlvals)-log(1-Correlvals));
%         
%     end
    
    MeanFisher{condition} = mean(Fishermatrix(:,:,:,condition),3);
end
%save('Fishermatrix_byload.mat','Fishermatrix')

pvector = [];

for i = 1:length(roi_names)
    for j = 1:length(roi_names)
        
        if i>j
            designmatrix = [];
            for subj = 1:length(subjects)
                for condition = 1:length(conditions)
                    designmatrix = [designmatrix ; [Fishermatrix(i,j,subj,condition) DATgenotype(subj) condition subj]];
                end
            end
            
            [SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(designmatrix,1);
            
            ANOVA_FVal_Load(i,j) = Fs{3};
            ANOVA_PVal_Load(i,j) = Ps{3};
            
            
            ANOVA_FVal_LoadxDAT(i,j) = Fs{4};
            ANOVA_PVal_LoadxDAT(i,j) = Ps{4};
   
            pvectorLoadxDAT = [pvector Ps{4}];
            pvector = [pvector Ps{3}];
   
            
            %                 for subj = 1:length(subjects)
            %                     texttowrite = [subjects{subj},'   ',runs{run},'   ',[roi_names{i} '_vs_' roi_names{j}],'   ',num2str(Fishermatrix{run}(i,j,subj))];
            %                     dlmwrite([analysis_dir outputfilename],texttowrite,'-append','delimiter','');
            %                 end
            
        elseif i == j
            MeanFisher{1}(i,j) = Inf; MeanFisher{2}(i,j) = Inf; MeanFisher{3}(i,j) = Inf; 
        else
            MeanFisher{1}(i,j) = 0; MeanFisher{2}(i,j) = 0; MeanFisher{3}(i,j) = 0; 
        end
        
    end
end



for condition = 1:length(conditions)
    imagesc(MeanFisher{condition},[-1,1])
    set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
    set(gca,'YTickLabel',roi_names);
    title([conditions{condition} ' Z Values']);
    
    %colormap autumn; autcmap = colormap; colormap winter; wintcmap = colormap;
    %combined = [flipdim(wintcmap,1); ones(1,3); autcmap];
    colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
    hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
    combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
    colormap(combined);
    %set(gca,'XGrid','on')
    %set(gca,'YGrid','on')
    colorbar
    
    set(gcf,'Nextplot','new');
end


alpha = .05;
%alpha = .05/231;
useFDR = 1;


%Load
%------------------------------------------------------------------------
[pID, pN] = FDR(pvector,alpha);
if isempty(pID); pID = 0; end

[pIDMax, pN] = FDR(pvector,.001);
if ~isempty(pIDMax); Flimit = finv(pIDMax,2,72); else Flimit = 7.61; end

FtestFDRBin = zeros(size(ANOVA_FVal_Load));
FtestFDRBin(find(ANOVA_PVal_Load<= pID)) = 1;


if useFDR; SigFtestVal = ANOVA_FVal_Load .* FtestFDRBin; else SigFtestVal = ANOVA_FVal_Load .* (ANOVA_PVal_Load <= alpha); end

imagesc(SigFtestVal,[0,Flimit])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title(['Load effect F Values']);

%colormap autumn; autcmap = colormap; colormap winter; wintcmap = colormap;
%combined = [flipdim(wintcmap,1); ones(1,3); autcmap];
colormap hot; 
% hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
% hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
% combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
% colormap(combined);
%set(gca,'XGrid','on')
%set(gca,'YGrid','on')
colorbar

set(gcf,'Nextplot','new');
%------------------------------------------------------------------------


%LoadxDAT
%------------------------------------------------------------------------
[pID, pN] = FDR(pvectorLoadxDAT,alpha);
if isempty(pID); pID = 0; end

[pIDMax, pN] = FDR(pvectorLoadxDAT,.001);
if ~isempty(pIDMax); Flimit = finv(pIDMax,2,72); else Flimit = 7.61; end

FtestFDRBin = zeros(size(ANOVA_FVal_LoadxDAT));
FtestFDRBin(find(ANOVA_PVal_LoadxDAT<= pID)) = 1;


if useFDR; SigFtestVal = ANOVA_FVal_LoadxDAT .* FtestFDRBin; else SigFtestVal = ANOVA_FVal_LoadxDAT .* (ANOVA_PVal_LoadxDAT <= alpha); end

imagesc(SigFtestVal,[0,Flimit])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title(['LoadxDAT effect F Values']);

%colormap autumn; autcmap = colormap; colormap winter; wintcmap = colormap;
%combined = [flipdim(wintcmap,1); ones(1,3); autcmap];
colormap hot; 
% hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
% hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
% combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
% colormap(combined);
%set(gca,'XGrid','on')
%set(gca,'YGrid','on')
colorbar

set(gcf,'Nextplot','new');
%------------------------------------------------------------------------






