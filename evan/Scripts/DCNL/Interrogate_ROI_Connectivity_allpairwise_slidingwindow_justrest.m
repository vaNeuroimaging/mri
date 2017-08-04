
warning off

directory = pwd;

basealpha = .05;
useBonf = 0;
useFDR = 0;

remove_FD = 1;
    FDthresh = .5;

windowsize = 30;

makenew = 0;

behaviorcorrel = 1;

outputfilename = 'SlidingConnectivityOutput.txt';

Datasavename = 'PairwiseData_constrained_withintask_slidingwindow_justrest';

analysis_dir = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Power_ROIs/';


addpath /fmri/data3/Evan/Gene-Rest-Nback/Scripts/LinStats2011/



subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%{'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
%{'269','110','189','199'};



genematrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Data/Genes.txt');

for sub = 1:length(subjects);
    COMTgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),3);
    DATgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),2);
    DRD4genotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),4);
    BDNFgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),5);
    SRTTgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),6);
end


behavname = {'RT','RT3B','Acc','Acc100%','Acc3B','ICV','ICV_3Back','Inattentive','H-I','BIS','STAI_Trait','Prescan_2B_Acc','Prescan_3B_Acc','Prescan_Overall_Acc'};

behavmatrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Behavior/AllMeasures.txt');

for behavnum = 1:length(behavname)
    for sub = 1:length(subjects)
        behavior{behavnum}(sub) = behavmatrix(find(behavmatrix(:,1)==str2num(subjects{sub})),behavnum+1);
    end
end


runs = {'FirstRest'};
    
%{'Rest','Nback'};

header = 'sw';

allrois = dir([analysis_dir '*_roi.mat']);
for i = 1:length(allrois)
roi_names{i} = allrois(i).name(1:end-8);
end


%roi_names = {'TNN-aDMN-vmPFC','TNN-aDMN-PCC','TNN-pDMN-Prec','TNN-pDMN-lAG','TNN-pDMN-rAG','TPN-Sal-dACC','TPN-Sal-raIns','TPN-Sal-laIns','TPN-rFPC-dlPFC','TPN-rFPC-IPL','TPN-lFPC-dlPFC','TPN-lFPC-IPL'};


if makenew ==0
    load([analysis_dir Datasavename]);
    
else
    delete([analysis_dir outputfilename]);
    fid = fopen([analysis_dir outputfilename],'at');
    fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Subject','Condition','ROI_Comparison','Fisher');
    fclose(fid);
    dlmwrite([analysis_dir outputfilename],' ','-append');
end
allpvector = [];


if makenew ==1
    
    for run = 1:length(runs)
        
        
        
        for subject = 1:length(subjects)
            
            clear timecourse residual_timecourse residuals Timecourses
            
            clear datamatrix;
            
            subjsid = subjects{subject};
            
            %data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/' runs{run} '/'];
            data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/' runs{run} '/'];
            
            motionparamfile = dir([data_dir 'rp*.txt']);
            motionparams = textread([data_dir motionparamfile.name]);
            
            clear roi_files;
            clear des_path;
            clear rois;
            clear des;
            clear mY;
            clear ttfile;
            clear imgfiles;
            clear n;
            clear P;
            clear raw_d;
            clear fmri_raw_data;
            clear numVox;
            
            imgfiles = dir([data_dir header '*.hdr']);
            m = size(imgfiles, 1);
            if ~isempty(imgfiles)
                for j=1:m
                    P(j, :) = [data_dir imgfiles(j).name];
                end
            else
                imgfiles = dir([data_dir header '*.nii']);
                m = size(imgfiles, 1);
                for j=1:m
                    P(j, :) = [data_dir imgfiles(j).name];
                end
            end
            
            
            
            if remove_FD == 1;
           
                FD = Calc_FD(data_dir);
           
                alltimepointstouse = find(FD<FDthresh);
           
                disp(['Subject ' subjsid ', ' runs{run} ': retained ' num2str(length(alltimepointstouse)) ' of ' num2str(size(P,1)) ' timepoints.'])
                
%                 if length(timepointstouse)>125;
%                     alltimepointstouse = timepointstouse(1:125);
%                 end
                
                P = P(alltimepointstouse,:);
                motionparams = motionparams(alltimepointstouse,:);
            
            end
            
            
            
                        
            disp([subjsid ' : ' runs{run} ' : CSF'])
            try
                LV_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_LV_roi.mat']);
            catch
                LV_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_CSF_roi.mat']);
            end
            [Y a b c] = getdata(LV_rois{1}, P,'l');
            LV_timecourse{run} = mean(Y,2);
            
            clear Y
            
            disp([subjsid ' : ' runs{run} ' : WM'])
            WM_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_WM_roi.mat']);
            [Y a b c] = getdata(WM_rois{1}, P,'l');
            WM_timecourse{run} = mean(Y,2);
            
            clear Y
            
            for roinum = 1:length(roi_names)
                disp([subjsid ' : ' runs{run} ' : ' roi_names{roinum}])
                roi = [analysis_dir, roi_names{roinum} , '_roi.mat'];
                %roi = [analysis_dir, roi_names{roinum}];
                rois = maroi('load_cell', roi);
                
                [Y a b c] = getdata(rois{1}, P,'l');
                
                timecourse{run,roinum} = mean(Y,2);
                
                for voxelnum = 1:size(Y,2)
                    
                    [b, bint, r] = regress(Y(:,voxelnum),[LV_timecourse{run} WM_timecourse{run} motionparams ones(length(LV_timecourse{run}),1)]);
                    
                    residuals(:,voxelnum) = r;
                    
                end
                
                residual_timepoint_means = mean(residuals,2);
                clear residuals filtered_residuals
                
                Timecourses(:,roinum) = residual_timepoint_means;
                
            end
            
            
            for windownum = 1:(length(Timecourses(:,1))-windowsize)
                correlvals = corrcoef(Timecourses(windownum:windownum+windowsize,:));
                Fisherwindowcorrelvectors{run}(windownum,:,:,subject) = .5*(log(1+correlvals)-log(1-correlvals));
            end
            Stdwindowcorrelvectors(:,:,subject,run) = std(Fisherwindowcorrelvectors{run}(:,:,:,subject));
            
  
            
        end
        
    end
    
    
end





clear FDRpvector

comparisoncounter = 0;


for i = 1:length(roi_names)
    for j = 1:length(roi_names)
        
        if i>j
            
            clear variablesinmodel modelmatrix m anovavals modelfit
            
            comparisoncounter = comparisoncounter + 1;
            
            
            if makenew ==1
                for run = 1:length(runs)
                    for subj = 1:length(subjects)
                        texttowrite = [subjects{subj},' ',runs{run},' ',[roi_names{i} '_vs_' roi_names{j}],' ',num2str(Stdwindowcorrelvectors(i,j,subj,run))];
                        dlmwrite([analysis_dir outputfilename],texttowrite,'-append','delimiter','');
                    end
                end
            end
            
%             [h,p,ci,stats] = ttest(Stdwindowcorrelvectors(i,j,:,2),Stdwindowcorrelvectors(i,j,:,1),basealpha);
%             RunTtestTVal(i,j) = stats.tstat;
%             RunTtestPVal(i,j) = p;
%             FDRpvector(comparisoncounter) = p;
            
            
            
            
        else
%             RunTtestTVal(i,j) = 0;
%             RunTtestPVal(i,j) = 1;
        end
        
    end
end
close all





if useFDR == 1;
    
    thisalpha = FDR(FDRpvector, basealpha);
    if isempty(thisalpha); thisalpha = 0; end
    
elseif useBonf == 1;
    
    thisalpha = basealpha / comparisoncounter;
    
else
    
    thisalpha = basealpha;
    
end

Tlimit = tinv((1-thisalpha/500),length(subjects)-2);
if Tlimit == inf
    Tlimit = 100000000000000000000;
end

% SigTtest = RunTtestPVal < thisalpha;
% SigTvals = RunTtestTVal .* SigTtest;
% 
% 
% imagesc(SigTvals,[-Tlimit,Tlimit])
% set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
% set(gca,'XTickLabel',roi_names);
% set(gca,'YTickLabel',roi_names);
% title([runs{2} ' > ' runs{1} ': T Values']);
% colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
% hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
% combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
% colormap(combined);
% colorbar
% 
% set(gcf,'Nextplot','new');
% 
% 
% 
% imagesc(RunTtestTVal,[-Tlimit,Tlimit])
% set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
% set(gca,'XTickLabel',roi_names);
% set(gca,'YTickLabel',roi_names);
% title([runs{2} ' > ' runs{1} ': Unthresholded T Values']);
% colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
% hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
% combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
% colormap(combined);
% colorbar
% 
% set(gcf,'Nextplot','new');



%correlations with behavior

if behaviorcorrel ==1
    for run=1:length(runs)
        for behavnum = 1:length(behavior)
                        
            CorrelBin{run} = zeros(length(roi_names));
            CorrelSig{run} = zeros(length(roi_names));
            CorrelVal{run} = zeros(length(roi_names));
            CorrelFDRBin{run} = zeros(length(roi_names));
            pvector = [];
            
            for i = 1:length(roi_names)
                for j = 1:length(roi_names)
                    
                    if i>j
                        
                        [R,P]=CORRCOEF(Stdwindowcorrelvectors(i,j,:,run),behavior{behavnum});
                        CorrelSig{run}(i,j) = P(2,1);
                        pvector(end+1) = P(2,1);
                        CorrelVal{run}(i,j) = R(2,1);
                        
                    else
                        CorrelSig{run}(i,j) = 0;
                    end
                    
                end
            end
            
            [pID, pN] = FDR(pvector,.05);
            if isempty(pID); pID = 0; end
            
            Rlimit = 1;
            
            CorrelFDRBin{run}(find(CorrelSig{run}>0 & CorrelSig{run}<= pID)) = 1;
            
            
            if useFDR; SigCorrelVal{run} = CorrelVal{run} .* CorrelFDRBin{run}; else; SigCorrelVal{run} = CorrelVal{run} .* (CorrelSig{run} < thisalpha) ; end
            
            imagesc(SigCorrelVal{run},[-Rlimit,Rlimit])
            set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
            set(gca,'XTickLabel',roi_names);
            set(gca,'YTickLabel',roi_names);
            title([runs{run} ' ' behavname{behavnum} ' R Values']);
            
            colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
            hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
            combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
            colormap(combined);
            colorbar
            
            set(gcf,'Nextplot','new');
            
        end
        
    end
end


% CorrelSig = zeros(length(roi_names));
% CorrelVal = zeros(length(roi_names));
%
% for i = 1:length(roi_names)
%  for j = 1:length(roi_names)
%
%  if i>j
%
%  [R,P]=CORRCOEF(Stdwindowcorrelvectors(i,j,:,1),Stdwindowcorrelvectors(i,j,:,2));
%  CorrelSig(i,j) = P(2,1);
%  CorrelVal(i,j) = R(2,1);
%
%  else
%  CorrelSig(i,j) = 0;
%  end
%
%  end
% end
%
%
% SigCorrelVal = CorrelVal .* (CorrelSig>0 & CorrelSig<= basealpha);
%
% Rlimit = 1;
%
% imagesc(SigCorrelVal,[-Rlimit,Rlimit])
% set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
% set(gca,'XTickLabel',roi_names);
% set(gca,'YTickLabel',roi_names);
% title(['Correlations between runs: R Values']);
%
% colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
% hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
% combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
% colormap(combined);
% colorbar
%
% set(gcf,'Nextplot','new');






clear valuestokeep
valuestokeep = zeros(length(roi_names));
for i = 1:length(roi_names)
    for j = 1:length(roi_names)
        if i>=j
            valuestokeep(i,j) = 1;
        end
    end
end

for run = 1:length(runs)
    
    imagesc(mean(Stdwindowcorrelvectors(:,:,:,run),3).*valuestokeep,[0,.3])
    set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
    set(gca,'YTickLabel',roi_names);
    title([runs{run} ' Z Values']);
    
    colormap hot;
    colorbar
    
    set(gcf,'Nextplot','new');
    
    
end


% MeanFisher_all = mean(mean(Stdwindowcorrelvectors,4),3);
% 
% imagesc(MeanFisher_all.*valuestokeep,[0,.3])
% set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
% set(gca,'XTickLabel',roi_names);
% set(gca,'YTickLabel',roi_names);
% title('Overall Z Values');
% 
% colormap hot;
% colorbar
% 
% set(gcf,'Nextplot','new');


if makenew ==1
    save([analysis_dir Datasavename],'Stdwindowcorrelvectors','Fisherwindowcorrelvectors')
end








