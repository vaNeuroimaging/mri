
warning off

directory = pwd;

basealpha = .05;
useBonf = 1;
useFDR = 0;
useMonteCarlo = 0;
simulations = 10000;

makenew = 0;

behaviorcorrel = 0;

outputfilename = 'ConnectivityOutput_filt_wholebrain.txt';
%'ConnectivityOutput_spm8.txt';

Datasavename = 'PairwiseData_filt_wholebrain';
%;'PairwiseData_constrained_withintask_noloadregression'

analysis_dir = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_77DAT_true.gica/ROIs/6mm/';
%'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_77DAT_true.gica/ROIs/6mm/';
%
%
% 
%

addpath /fmri/data3/Evan/Gene-Rest-Nback/Scripts/LinStats2011/



subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
%,'181','182'


genematrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Data/Genes.txt');

for sub = 1:length(subjects);
    COMTgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),3);
    DATgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),2);
    DRD4genotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),4);
    BDNFgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),5);
    SRTTgenotype(sub) = genematrix(find(genematrix==str2num(subjects{sub})),6);
end

behavname = {'RT','RT3B','Acc','Acc100%','Acc3B','ICV','ICV_3Back','Inattentive','H-I','BIS'};

behavmatrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Behavior/AllMeasures.txt');

for behavnum = 1:length(behavname)
    for sub = 1:length(subjects)
        behavior{behavnum}(sub) = behavmatrix(find(behavmatrix(:,1)==str2num(subjects{sub})),behavnum+1);
    end
end


runs = {'Rest','Nback'};
%{'Rest','Nback'};
% conditions = {'2Back','3Back'};
% condtimes = {[5 115 159] [27 71 181]};
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


header = 'fsw';


allrois = dir([analysis_dir '*_roi.mat']);
for i = 1:length(allrois)
    roi_names{i} = allrois(i).name(1:end-8);
end
roi_names = {'rFPC-rPar'    'rFPC-rdlPFC'    'rFPC-rvlPFC' 'lFPC-ldlPFC'   'DA-rPar' 'DA-lPar' 'Sal-raIns' 'Sal-laIns'  'Sal-ramfg'  'Sal-lamfg'  'Striatum-r' 'Striatum-l'  'aDMN-vmPFC' 'aDMN-PCC' 'pDMN-rAG' 'pDMN-lAG'};
%'bilatfront-l','bilatfront-r',
%{'rFPC-rdlPFC','rFPC-rPar','lFPC-ldlPFC','lFPC-lPar','Sal-raIns','Sal-laIns','Sal-ramfg','Sal-lamfg','DA-rPar','DA-lPar','Striatum-r','Striatum-l','aDMN-vmPFC','pDMN-Prec','pDMN-rAG','pDMN-lAG'};
%'rFPC-rvlPFC','lFPC-SMA',
%{'rFPC-rdlPFC'    'rFPC-rvlPFC'    'rFPC-rPar'    'lFPC-ldlPFC'  'lFPC-lPar'  'Sal-raIns'  'Sal-laIns'    'Sal-raMFG'   'Sal-laMFG'   'aDMN-vmPFC'  'aDMN-PCC'  'aDMN-lAG'  'pDMN-rPrec'  'pDMN-lPrec'  'pDMN-lAG' };



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
            
            clear timecourse residual_timecourse residuals
            
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
            if isempty(imgfiles)
                imgfiles = dir([data_dir header '*.nii']);
            end
            m = size(imgfiles, 1);
            for j=1:m
                P(j, :) = [data_dir imgfiles(j).name];
            end
            
            if strcmp(runs{run},'Nback')
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
            
            
            disp([subjsid ' : ' runs{run} ' : Wholebrain'])
            WB_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_roi.mat']);
            [Y a b c] = getdata(WB_rois{1}, P,'l');
            WB_timecourse{run} = mean(Y,2);
            
            clear Y
            
            for roinum = 1:length(roi_names)
                disp([subjsid ' : ' runs{run} ' : ' roi_names{roinum}])
                roi = [analysis_dir, roi_names{roinum} , '_roi.mat'];
                %roi = [analysis_dir, roi_names{roinum}];
                rois = maroi('load_cell', roi);
                
                [Y a b c] = getdata(rois{1}, P,'l');
                
                timecourse{run,roinum} = mean(Y,2);
                
                for voxelnum = 1:size(Y,2)
                    %[b, bint, r] = regress(Y(:,voxelnum),[LV_timecourse{run} WM_timecourse{run} motionparams ones(length(LV_timecourse{run}),1)]);
                    if strcmp(runs{run},'Nback')
                        [b, bint, r] = regress(Y(:,voxelnum),[LV_timecourse{run} WM_timecourse{run} WB_timecourse{run} motionparams taskregressor ones(length(LV_timecourse{run}),1)]);
                    else
                        [b, bint, r] = regress(Y(:,voxelnum),[LV_timecourse{run} WM_timecourse{run} WB_timecourse{run} motionparams ones(length(LV_timecourse{run}),1)]);
                    end
                    
                    
                    residuals(:,voxelnum) = r;
                    
                end
                
                residual_timepoint_means = mean(residuals,2);
                clear residuals filtered_residuals
                
                Timecourses{run}(:,roinum) = residual_timepoint_means;
                
            end
            
            
            Correlvals = corrcoef(Timecourses{run}(:,:));
            Fishermatrix(:,:,subject,run) = .5*(log(1+Correlvals)-log(1-Correlvals));
            
        end
        
    end
    
    
    FDRpvectorRun = [];
    FDRpvectorDAT = [];
    FDRpvectorRunxDAT = [];
    FDRpvectorModel = [];
    
end






% if useMonteCarlo==1
%     for i = 1:length(roi_names)
%         for j = 1:length(roi_names)
%             if i>j
%                 mu = mean(mean(Fishermatrix(i,j,:,:)));
%                 sigma = std(reshape(Fishermatrix(i,j,:,:),[1,size(Fishermatrix,3)*size(Fishermatrix,4)]));
%
%                 SimulatedData(i,j,:,:,:) = normrnd(mu,sigma,size(Fishermatrix,3),size(Fishermatrix,4),simulations);
%
%             end
%         end
%     end
% end




pvector_MonteCarlo_Run = [];
pvector_MonteCarlo_DAT = [];
pvector_MonteCarlo_RunxDAT = [];


comparisoncounter = 0;

if useMonteCarlo==1
    h = waitbar(0,'Monte-Carlo Simulating');
end

for i = 1:length(roi_names)
    for j = 1:length(roi_names)
        
        if i>j
            
            clear variablesinmodel modelmatrix m anovavals modelfit
            
            comparisoncounter = comparisoncounter + 1;
            for run = 1:length(runs)
                [h,p,ci,stats] = ttest2(Fishermatrix(i,j,find(DATgenotype==2),run),Fishermatrix(i,j,find(DATgenotype==1),run),basealpha);
                DATTtestVal_indRun{run}(i,j) = stats.tstat;
            end
            
            if makenew==1

                for run = 1:length(runs)
                    for subj = 1:length(subjects)
                        texttowrite = [subjects{subj},'   ',runs{run},'   ',[roi_names{i} '_vs_' roi_names{j}],'   ',num2str(Fishermatrix(i,j,subj,run))];
                        dlmwrite([analysis_dir outputfilename],texttowrite,'-append','delimiter','');
                    end
                end
                
                [h,p,ci,stats] = ttest(Fishermatrix(i,j,:,2),Fishermatrix(i,j,:,1),basealpha);
                RunTtestVal(i,j) = stats.tstat;
                
                [h,p,ci,stats] = ttest2(reshape(Fishermatrix(i,j,find(DATgenotype==2),:),[1,length(find(DATgenotype==2))*2]),reshape(Fishermatrix(i,j,find(DATgenotype==1),:),[1,length(find(DATgenotype==1))*2]),basealpha);
                DATTtestVal(i,j) = stats.tstat;
                
                
                
                
                designmatrix = [];
                for subj = 1:length(subjects)
                    for run = 1:length(runs)
                        designmatrix = [designmatrix ; [Fishermatrix(i,j,subj,run) DATgenotype(subj) run subj]];
                    end
                end
                
                variablesinmodel = Vars(designmatrix(:,2),designmatrix(:,3),designmatrix(:,4),'type',[1 1 3]);
                modelmatrix = [0 0 0; eye(3); 1 1 0];
                m = Model(designmatrix(:,1),variablesinmodel,modelmatrix);
                
                anovavals = m.anova;
                
                ANOVA_FVal_DAT(i,j) = anovavals.abs(2,1);
                ANOVA_PVal_DAT(i,j) = anovavals.abs(2,4);
                FDRpvectorDAT = [FDRpvectorDAT anovavals.abs(2,4)];
                
                ANOVA_FVal_Run(i,j) = anovavals.abs(3,1);
                ANOVA_PVal_Run(i,j) = anovavals.abs(3,4);
                FDRpvectorRun = [FDRpvectorRun anovavals.abs(3,4)];
                
                ANOVA_FVal_RunxDAT(i,j) = anovavals.abs(4,1);
                ANOVA_PVal_RunxDAT(i,j) = anovavals.abs(4,4);
                FDRpvectorRunxDAT = [FDRpvectorRunxDAT anovavals.abs(4,4)];
                
                modelfit = m.fit;
                
                ANOVA_FVal_Model(i,j) = modelfit.abs(1,1);
                ANOVA_PVal_Model(i,j) =  modelfit.abs(1,4);
                FDRpvectorModel = [FDRpvectorModel modelfit.abs(1,4)];
                
                
                
                %             [SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(designmatrix,1);
                %
                %             ANOVA_FVal_DAT(i,j) = Fs{1};
                %             ANOVA_PVal_DAT(i,j) = Ps{1};
                %             FDRpvectorDAT = [FDRpvectorDAT Ps{1}];
                %
                %             ANOVA_FVal_Run(i,j) = Fs{3};
                %             ANOVA_PVal_Run(i,j) = Ps{3};
                %             FDRpvectorRun = [FDRpvectorRun Ps{3}];
                %
                %             ANOVA_FVal_RunxDAT(i,j) = Fs{4};
                %             ANOVA_PVal_RunxDAT(i,j) = Ps{4};
                %             FDRpvectorRunxDAT = [FDRpvectorRunxDAT Ps{4}];
                
                %             regressiondesignmatrix = [designmatrix designmatrix(:,2).*designmatrix(:,3) ones(size(designmatrix,1),1)];
                %             [B,BINT,R,RINT,STATS] = REGRESS(regressiondesignmatrix(:,1),regressiondesignmatrix(:,2:end));
                %             Regression_Pval(i,j) = STATS(3);
                %             Regression_Fval(i,j) = STATS(2);
                %             regressionpvector = [regressionpvector STATS(3)];
                
                if useMonteCarlo==1
                    for simnum = 1:simulations
                        waitbar(((comparisoncounter-1)*simulations + simnum) / ((length(roi_names)^2 - length(roi_names))*simulations/2));
                        designmatrix = [];
                        subjperm = randperm(length(subjects));
                        for subj = 1:length(subjects)
                            runperm = randperm(2);
                            for run = 1:length(runs)
                                designmatrix = [designmatrix ; [Fishermatrix(i,j,subj,run) DATgenotype(subjperm(subj)) runperm(run) subjperm(subj)]];
                                %designmatrix = [designmatrix ; [SimulatedData(i,j,subj,run,simnum) DATgenotype(subj) run subj]];
                            end
                        end
                        [SSQs, DFs, MSQs, Fs, Ps]=mixed_between_within_anova(designmatrix,1);
                        pvector_MonteCarlo_Run = [pvector_MonteCarlo_Run Ps{3}];
                        pvector_MonteCarlo_DAT = [pvector_MonteCarlo_DAT Ps{1}];
                        pvector_MonteCarlo_RunxDAT = [pvector_MonteCarlo_RunxDAT Ps{4}];
                    end
                end
            end
            
        else
            ANOVA_FVal_DAT(i,j) = 0;
            ANOVA_PVal_DAT(i,j) = 0;
            ANOVA_FVal_RunxDAT(i,j) = 0;
            ANOVA_PVal_RunxDAT(i,j) = 0;
            ANOVA_FVal_Run(i,j) = 0;
            ANOVA_PVal_Run(i,j) = 0;
            RunTtestVal(i,j) = 0;
            DATTtestVal(i,j) = 0;
            ANOVA_FVal_Model(i,j) = 0;
            ANOVA_PVal_Model(i,j) =  0;
            DATTtestVal_indRun{1}(i,j) = 0;
            DATTtestVal_indRun{2}(i,j) = 0;
        end
        
    end
end
close all



%Testing effects post-hoc, after overall model significance is established
%%%%% Model Significance %%%%%%
%------------------------------------------------------------------------



if useFDR == 1;
    
    thisalpha = FDR(FDRpvectorModel, basealpha);
    if isempty(thisalpha); thisalpha = 0; end
    
elseif useBonf == 1;
    
    thisalpha = basealpha / comparisoncounter;
    
    % elseif useMonteCarlo == 1;
    %     pvector_MonteCarlo_Run = sort(pvector_MonteCarlo_Run);
    %     thisalpha = pvector_MonteCarlo_Run(comparisoncounter*simulations*basealpha)/comparisoncounter
    
else
    
    thisalpha = basealpha;
    
end

Flimit = Finv(1-(thisalpha/50),1,length(subjects)-2);

SigANOVA_PVal_Model = ANOVA_PVal_Model < thisalpha;

SigANOVA_FVal_Model = ANOVA_FVal_Model .* SigANOVA_PVal_Model;

% imagesc(SigANOVA_FVal_Model,[0,Flimit])
% set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
% set(gca,'XTickLabel',roi_names);
% set(gca,'YTickLabel',roi_names);
% title(['Overall Model significance']);
% 
% colormap hot;
% colorbar
% set(gcf,'Nextplot','new');



%%%%% STATExDAT EFFECT %%%%%%
%------------------------------------------------------------------------

Flimit = Finv(1-(basealpha/50),1,length(subjects)-2);

SigANOVA_PVal_RunxDAT = (ANOVA_PVal_RunxDAT < basealpha) .* SigANOVA_PVal_Model;

RunxDATSigFtestVal = ANOVA_FVal_RunxDAT .* SigANOVA_PVal_RunxDAT;

imagesc(RunxDATSigFtestVal,[0,Flimit])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title(['RunxDAT effect F Values']);

colormap hot;
colorbar
set(gcf,'Nextplot','new');


%%%%% STATE EFFECT %%%%%%
%------------------------------------------------------------------------
Tlimit = -tinv(1-(basealpha/500000),length(subjects)-2);

SigANOVA_PVal_Run = ((ANOVA_PVal_Run < basealpha) .* SigANOVA_PVal_Model);%.* (~SigANOVA_PVal_RunxDAT);

RunSigTtestVal = RunTtestVal .* SigANOVA_PVal_Run;

imagesc(RunSigTtestVal,[Tlimit,-Tlimit])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title([runs{2} ' vs ' runs{1} ': T Values']);
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
colorbar

set(gcf,'Nextplot','new');


%%%%% DAT EFFECT %%%%%%
%------------------------------------------------------------------------


SigANOVA_PVal_DAT = ((ANOVA_PVal_DAT < basealpha) .* SigANOVA_PVal_Model); %.* (~SigANOVA_PVal_RunxDAT);

DATSigTtestVal = DATTtestVal .* SigANOVA_PVal_DAT;

imagesc(DATSigTtestVal,[Tlimit,-Tlimit])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title(['DAT Main effect T Values']);
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
colorbar

set(gcf,'Nextplot','new');





%correlations with behavior

if behaviorcorrel ==1
    
    for behavnum = 1:length(behavior)
        for run=1:length(runs)
            
            CorrelBin{run} = zeros(length(roi_names));
            CorrelSig{run} = zeros(length(roi_names));
            CorrelVal{run} = zeros(length(roi_names));
            CorrelFDRBin{run} = zeros(length(roi_names));
            pvector = [];
            
            for i = 1:length(roi_names)
                for j = 1:length(roi_names)
                    
                    if i>j
                        
                        [R,P]=CORRCOEF(Fishermatrix(i,j,:,run),behavior{behavnum});
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
        
        CorrelBinRundiff = zeros(length(roi_names));
        CorrelSigRundiff = zeros(length(roi_names));
        CorrelValRundiff = zeros(length(roi_names));
        CorrelFDRBinRundiff = zeros(length(roi_names));
        pvector = [];
        
        for i = 1:length(roi_names)
            for j = 1:length(roi_names)
                
                if i>j
                    
                    [R,P]=CORRCOEF((Fishermatrix(i,j,:,2) - Fishermatrix(i,j,:,1)),behavior{behavnum});
                    CorrelSigRundiff(i,j) = P(2,1);
                    pvector(end+1) = P(2,1);
                    CorrelValRundiff(i,j) = R(2,1);
                    
                else
                    CorrelSigRundiff(i,j) = 0;
                end
                
            end
        end
        
        [pID, pN] = FDR(pvector,.05);
        if isempty(pID); pID = 0; end
        
        Rlimit = 1;
        
        CorrelFDRBinRundiff(find(CorrelSigRundiff>0 & CorrelSigRundiff<= pID)) = 1;
        
        
        if useFDR; SigCorrelValRundiff = CorrelValRundiff .* CorrelFDRBinRundiff; else; SigCorrelValRundiff = CorrelValRundiff .* (CorrelSigRundiff < thisalpha) ; end
        
        imagesc(SigCorrelValRundiff,[-Rlimit,Rlimit])
        set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
        set(gca,'XTickLabel',roi_names);
        set(gca,'YTickLabel',roi_names);
        title(['Nback-Rest ' behavname{behavnum} ' R Values']);
        
        colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
        hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
        combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
        colormap(combined);
        colorbar
        
        set(gcf,'Nextplot','new');
        
        
    end
end


CorrelSig = zeros(length(roi_names));
CorrelVal = zeros(length(roi_names));

for i = 1:length(roi_names)
            for j = 1:length(roi_names)
                
                if i>j
                    
                    [R,P]=CORRCOEF(Fishermatrix(i,j,:,1),Fishermatrix(i,j,:,2));
                    CorrelSig(i,j) = P(2,1);
                    CorrelVal(i,j) = R(2,1);
                    
                else
                    CorrelSig(i,j) = 0;
                end
                
            end
end



%Correlation between runs

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
    
    MeanFisher_910{run} = mean(Fishermatrix(:,:,find(DATgenotype==1),run),3);
    MeanFisher_1010{run} =mean(Fishermatrix(:,:,find(DATgenotype==2),run),3);
    
    
    imagesc(MeanFisher_910{run}.*valuestokeep,[-1,1])
    set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
    set(gca,'YTickLabel',roi_names);
    title([runs{run} ' 9/10 Z Values']);
    
    colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
    hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
    combined = [flipdim(coolmap,1); hotmap];
    colormap(combined);
    colorbar
    
    set(gcf,'Nextplot','new');
    
    
    imagesc(MeanFisher_1010{run}.*valuestokeep,[-1,1])
    set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
    set(gca,'YTickLabel',roi_names);
    title([runs{run} ' 10/10 Z Values']);
    
    colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
    hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
    combined = [flipdim(coolmap,1); hotmap];
    colormap(combined);
    colorbar
    
    set(gcf,'Nextplot','new');
    
    
    
    
    %Ttest value
%     imagesc(DATTtestVal_indRun{run}.*valuestokeep,[Tlimit,-Tlimit])
%     set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
%     set(gca,'XTickLabel',roi_names);
%     set(gca,'YTickLabel',roi_names);
%     title([runs{run} ' DAT T-test T Values']);
%     
%     colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
%     hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
%     combined = [flipdim(coolmap,1); hotmap];
%     colormap(combined);
%     colorbar
%     
%     set(gcf,'Nextplot','new');
    
    
    
    
    
    
end


MeanFisher_all = mean(mean(Fishermatrix,4),3);

imagesc(MeanFisher_all.*valuestokeep,[-1,1])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title('Overall Z Values');

colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); hotmap];
colormap(combined);
colorbar

set(gcf,'Nextplot','new');


%
%
%
%
%
%     [pID, pN] = FDR(allpvector,.05);
%     if isempty(pID); pID = 0; end
%
%     [pIDMax, pN] = FDR(allpvector,.001);
%     if ~isempty(pIDMax) & useFDR ==1; Tlimit = -tinv(1-pIDMax,length(subjects)-2); elseif useFDR == 1 & pID>0; Tlimit = -tinv(1-pID/20,length(subjects)-2); else Tlimit = -tinv(1-alpha/20,length(subjects)-2); end
%
%     TtestFDRBin{run}(find(TtestSig{run}>0 & TtestSig{run}<= pID)) = 1;
%
%
%     if useFDR; SigTtestVal{run} = TtestVal{run} .* TtestFDRBin{run}; else; SigTtestVal{run} = TtestVal{run} .* TtestBin{run}; end
%
%     imagesc(SigTtestVal{run},[Tlimit,-Tlimit])
%     set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
%     set(gca,'XTickLabel',roi_names);
%     set(gca,'YTickLabel',roi_names);
%     title([runs{run} ' DAT T Values']);
%
%     %colormap autumn; autcmap = colormap; colormap winter; wintcmap = colormap;
%     %combined = [flipdim(wintcmap,1); ones(1,3); autcmap];
%     colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
%     hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
%     combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
%     colormap(combined);
%     %set(gca,'XGrid','on')
%     %set(gca,'YGrid','on')
%     colorbar
%
%     set(gcf,'Nextplot','new');
%
%
%
%
%     imagesc(TtestVal{run},[Tlimit,-Tlimit])
%     set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
%     set(gca,'XTickLabel',roi_names);
%     set(gca,'YTickLabel',roi_names);
%     title([runs{run} ' DAT unthresholded T Values']);
%
%     %colormap autumn; autcmap = colormap; colormap winter; wintcmap = colormap;
%     %combined = [flipdim(wintcmap,1); ones(1,3); autcmap];
%     colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
%     hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
%     %combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
%     combined = [flipdim(coolmap,1); hotmap];
%     colormap(combined);
%     %set(gca,'XGrid','on')
%     %set(gca,'YGrid','on')
%     colorbar
%
%     set(gcf,'Nextplot','new');
%
%
% end
%
% RunTtestBin = zeros(length(roi_names));
% RunTtestSig = zeros(length(roi_names));
% RunTtestVal = zeros(length(roi_names));
% RunTtestFDRBin = zeros(length(roi_names));


if makenew ==1
    save([analysis_dir Datasavename],'Fishermatrix','ANOVA_FVal_DAT','ANOVA_FVal_Model','ANOVA_FVal_Run','ANOVA_FVal_RunxDAT','ANOVA_PVal_DAT','ANOVA_PVal_Model','ANOVA_PVal_Run','ANOVA_PVal_RunxDAT','FDRpvectorDAT','FDRpvectorModel','FDRpvectorRun','FDRpvectorRunxDAT','RunTtestVal','DATTtestVal')
    %    save([analysis_dir 'PairwiseData_taskregressor'],'Fishermatrix');
end








