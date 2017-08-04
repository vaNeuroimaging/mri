
warning off

directory = pwd;

basealpha = .05;
useBonf = 1;

makenew = 0;

behaviorcorrel = 0;

outputfilename = 'ConnectivityOutput_filt.txt';
%This contains values that can be imported into e.g. excel to make graphs
%and such

Datasavename = 'PairwiseData_filt';
%This is a data file which gets saved so you don't have to rerun all the
%connectivity if you wnat to see the results again


analysis_dir = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_61DAT_filt.gica/ROIs/6mm/';
%Where the ROIs are saved, and where the outputs will be spit

addpath /fmri/data3/Evan/Gene-Rest-Nback/Scripts/LinStats2011/



subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327'};


%you need a "Group" variable that lists the group ID of each subject.  
Group = [1 1 2 1 2 2 1 2 1]; %etc.




%For correlating behavior with connectivity
behavname = {'RT','RT3B','Acc','Acc100%','Acc3B','ICV','ICV_3Back','Inattentive','H-I','BIS'};

behavmatrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Behavior/AllMeasures.txt');

for behavnum = 1:length(behavname)
    for sub = 1:length(subjects)
        behavior{behavnum}(sub) = behavmatrix(find(behavmatrix(:,1)==str2num(subjects{sub})),behavnum+1);
    end
end


runs = {'Rest','Nback'};


header = 'fsw';


%This just loads all rois in the analysis directory
allrois = dir([analysis_dir '*_roi.mat']);
for i = 1:length(allrois)
    roi_names{i} = allrois(i).name(1:end-8);
end
%roi_names = {'rFPC-rdlPFC' 'rFPC-rvlPFC' 'rFPC-rPar' 'sal-SMA' 'sal-laIns' 'sal-raIns' 'sal-lmfg' 'sal-rmfg' 'str-l' 'str-r' 'DMN-vmPFC' 'DMN-Prec' 'DMN-lAG'};



%If the makenew variable = 0, then it won't redo all the connectivities,
%but just show you the results again.  This is because it ain't quick to
%rerun.
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
            
            %Location of data for this subject
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
                        
            disp([subjsid ' : ' runs{run} ' : CSF'])
            
            %Load and extract data from CSF ROI
            LV_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_CSF_roi.mat']);
            [Y a b c] = getdata(LV_rois{1}, P,'l');
            LV_timecourse{run} = mean(Y,2);
            
            clear Y
            
            %Load and extract data from WM ROI
            disp([subjsid ' : ' runs{run} ' : WM'])
            WM_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_WM_roi.mat']);
            [Y a b c] = getdata(WM_rois{1}, P,'l');
            WM_timecourse{run} = mean(Y,2);
            
            clear Y
            
            for roinum = 1:length(roi_names)
                disp([subjsid ' : ' runs{run} ' : ' roi_names{roinum}])
                
                %Load and extract data from this ROI
                roi = [analysis_dir, roi_names{roinum} , '_roi.mat'];
                rois = maroi('load_cell', roi);
                [Y a b c] = getdata(rois{1}, P,'l');
                Timecourses{run}(:,roinum) = mean(Y,2);
                
            end
            
            %Run correlations and convert to Fisher Z values
            Correlvals = partialcorr(Timecourses{run},[LV_timecourse{run} WM_timecourse{run} motionparams]);
            Fishermatrix(:,:,subject,run) = .5*(log(1+Correlvals)-log(1-Correlvals));
            
        end
        
    end
    

    
end




comparisoncounter = 0;


%Run statistical tests

for i = 1:length(roi_names)
    for j = 1:length(roi_names)
        
        if i>j
            
            clear variablesinmodel modelmatrix m anovavals modelfit
            
            %Group T tests for each run
            comparisoncounter = comparisoncounter + 1;
            for run = 1:length(runs)
                [h,p,ci,stats] = ttest2(Fishermatrix(i,j,find(Group==2),run),Fishermatrix(i,j,find(Group==1),run),basealpha);
                GroupTtestVal_indRun{run}(i,j) = stats.tstat;
            end
            
            if makenew==1

                
                
                for run = 1:length(runs)
                    for subj = 1:length(subjects)
                        texttowrite = [subjects{subj},'   ',runs{run},'   ',[roi_names{i} '_vs_' roi_names{j}],'   ',num2str(Fishermatrix(i,j,subj,run))];
                        dlmwrite([analysis_dir outputfilename],texttowrite,'-append','delimiter','');
                    end
                end
                
                %paired t-tests looking for effect of run
                [h,p,ci,stats] = ttest(Fishermatrix(i,j,:,2),Fishermatrix(i,j,:,1),basealpha);
                RunTtestVal(i,j) = stats.tstat;
                
                %2-sample t-tests looking for effect of group, across runs
                [h,p,ci,stats] = ttest2(reshape(Fishermatrix(i,j,find(Group==2),:),[1,length(find(Group==2))*2]),reshape(Fishermatrix(i,j,find(Group==1),:),[1,length(find(Group==1))*2]),basealpha);
                GroupTtestVal(i,j) = stats.tstat;
                
                
                
                
                %Run x group anova
                
                designmatrix = [];
                for subj = 1:length(subjects)
                    for run = 1:length(runs)
                        designmatrix = [designmatrix ; [Fishermatrix(i,j,subj,run) Group(subj) run subj]];
                    end
                end
                
                variablesinmodel = Vars(designmatrix(:,2),designmatrix(:,3),designmatrix(:,4),'type',[1 1 3]);
                modelmatrix = [0 0 0; eye(3); 1 1 0];
                m = Model(designmatrix(:,1),variablesinmodel,modelmatrix);
                
                anovavals = m.anova;
                
                ANOVA_FVal_Group(i,j) = anovavals.abs(2,1);
                ANOVA_PVal_Group(i,j) = anovavals.abs(2,4);
                
                ANOVA_FVal_Run(i,j) = anovavals.abs(3,1);
                ANOVA_PVal_Run(i,j) = anovavals.abs(3,4);
                
                ANOVA_FVal_RunxGroup(i,j) = anovavals.abs(4,1);
                ANOVA_PVal_RunxGroup(i,j) = anovavals.abs(4,4);
                
                modelfit = m.fit;
                
                ANOVA_FVal_Model(i,j) = modelfit.abs(1,1);
                ANOVA_PVal_Model(i,j) =  modelfit.abs(1,4);
                
                
                
                
            end
            
        else
            ANOVA_FVal_Group(i,j) = 0;
            ANOVA_PVal_Group(i,j) = 0;
            ANOVA_FVal_RunxGroup(i,j) = 0;
            ANOVA_PVal_RunxGroup(i,j) = 0;
            ANOVA_FVal_Run(i,j) = 0;
            ANOVA_PVal_Run(i,j) = 0;
            RunTtestVal(i,j) = 0;
            GroupTtestVal(i,j) = 0;
            ANOVA_FVal_Model(i,j) = 0;
            ANOVA_PVal_Model(i,j) =  0;
            GroupTtestVal_indRun{1}(i,j) = 0;
            GroupTtestVal_indRun{2}(i,j) = 0;
        end
        
    end
end
close all




%Get the p value to use
if useBonf == 1;
    
    thisalpha = basealpha / comparisoncounter;
    
else
    
    thisalpha = basealpha;
    
end





%%%%% STATExGroup EFFECT %%%%%%
%------------------------------------------------------------------------

Flimit = Finv(1-(thisalpha/50),1,length(subjects)-2);

SigANOVA_PVal_RunxGroup = (ANOVA_PVal_RunxGroup < thisealpha);

RunxGroupSigFtestVal = ANOVA_FVal_RunxGroup .* SigANOVA_PVal_RunxGroup;


%Display it

imagesc(RunxGroupSigFtestVal,[0,Flimit])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title(['RunxGroup effect F Values']);

colormap hot;
colorbar
set(gcf,'Nextplot','new');



%%%%% STATE EFFECT %%%%%%
%------------------------------------------------------------------------
Tlimit = -tinv(1-(thisalpha/50),length(subjects)-2);

SigANOVA_PVal_Run = (ANOVA_PVal_Run < thisalpha);

RunSigTtestVal = RunTtestVal .* SigANOVA_PVal_Run;

%Display it

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


%%%%% Group EFFECT %%%%%%
%------------------------------------------------------------------------


SigANOVA_PVal_Group = (ANOVA_PVal_Group < thisalpha);

GroupSigTtestVal = GroupTtestVal .* SigANOVA_PVal_Group;

%Display it

imagesc(GroupSigTtestVal,[Tlimit,-Tlimit])
set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
set(gca,'XTickLabel',roi_names);
set(gca,'YTickLabel',roi_names);
title(['Group Main effect T Values']);
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
            
            if isempty(pID); pID = 0; end
            
            Rlimit = 1;
            
            
            
            
            SigCorrelVal{run} = CorrelVal{run} .* (CorrelSig{run} < thisalpha);
            
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
        
        Rlimit = 1;
        
        SigCorrelValRundiff = CorrelValRundiff .* (CorrelSigRundiff < thisalpha);
        
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






%Display average correlation matrices for each group in each run

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
    
    MeanFisher_Group1{run} = mean(Fishermatrix(:,:,find(Group==1),run),3);
    MeanFisher_Group2{run} =mean(Fishermatrix(:,:,find(Group==2),run),3);
    
    
    imagesc(MeanFisher_Group1{run}.*valuestokeep,[-1,1])
    set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
    set(gca,'YTickLabel',roi_names);
    title([runs{run} ' Group1 Z Values']);
    
    colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
    hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
    combined = [flipdim(coolmap,1); hotmap];
    colormap(combined);
    colorbar
    
    set(gcf,'Nextplot','new');
    
    
    imagesc(MeanFisher_Group2{run}.*valuestokeep,[-1,1])
    set(gca,'XTick',[1:length(roi_names)]); set(gca,'YTick',[1:length(roi_names)]);
    set(gca,'XTickLabel',roi_names);
    set(gca,'YTickLabel',roi_names);
    title([runs{run} ' Group2 Z Values']);
    
    colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
    hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
    combined = [flipdim(coolmap,1); hotmap];
    colormap(combined);
    colorbar
    
    set(gcf,'Nextplot','new');
    
    
    
    
    
    
    
    
    
end


%And display overall average Z values across all groups and runs

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



if makenew ==1
    save([analysis_dir Datasavename],'Fishermatrix','ANOVA_FVal_Group','ANOVA_FVal_Model','ANOVA_FVal_Run','ANOVA_FVal_RunxGroup','ANOVA_PVal_Group','ANOVA_PVal_Model','ANOVA_PVal_Run','ANOVA_PVal_RunxGroup','RunTtestVal','GroupTtestVal')
    %    save([analysis_dir 'PairwiseData_taskregressor'],'Fishermatrix');
end








