MSCnames = {'MSC01','MSC02','MSC03','MSC04'};
ncortverts = 59412;

sizethresh = 300;

names = {'DMN','Vis','FP','','DA','','VA','Sal','CO','MH','MM','Aud','','','PERN','White'};
oknetworknums = [1 2 3 5 7 8 9 10 11 12 15 16];

%PerceptronData = cifti_read('/data/cn4/evan/RestingState/Ind_variability/MSC/Perceptron_kden_results_MSC_clean50mm.dtseries.nii');
InfomapData = cifti_read('/data/cn4/evan/RestingState/Ind_variability/MSC/Consensus_Infomap_MSC.dtseries.nii');
TemplateData = cifti_read('/data/cn4/evan/RestingState/Ind_variability/MSC/Templatematch_dice_bysubject_betterkden.dtseries.nii');
Infomapinds = zeros(size(InfomapData));
for num = oknetworknums
    Infomapinds = Infomapinds + (InfomapData==num);
end
InfomapData = InfomapData .* Infomapinds;
InfomapData(isnan(InfomapData)) = 0;

%PminIcorrelmaps = zeros(66697,0);
%PminTcorrelmaps = zeros(66697,0);
IminTcorrelmaps = zeros(66697,0);

%PminIlabels = [];
%PminTlabels = [];
IminTlabels = [];

for s = 1:length(MSCnames)
    disp(['Subject ' num2str(s)])
    SAnames = {[MSCnames{s} '.L.midthickness.32k_fs_LR_surfaceareas.func.gii'],[MSCnames{s} '.R.midthickness.32k_fs_LR_surfaceareas.func.gii']};
%     PminI = (PerceptronData(:,s) - InfomapData(:,s)) .* (PerceptronData(:,s) > 0) .* (InfomapData(:,s) > 0);
%     IDs = unique(PminI); IDs(IDs==0) = [];
%     PminIclustereddiffs = zeros(size(PminI,1),0);
%     for ID = IDs(:)'
%          PminIclustereddiffs = [PminIclustereddiffs metric_cluster_cifti_surfacearea(PminI,ID-.001,ID+.001,sizethresh,SAnames) ];
%     end
%     
%     PminT = (PerceptronData(:,s) - TemplateData(:,s)) .* (PerceptronData(:,s) > 0) .* (TemplateData(:,s) > 0);
%     IDs = unique(PminT); IDs(IDs==0) = [];
%     PminTclustereddiffs = zeros(size(PminT,1),0);
%     for ID = IDs(:)'
%          PminTclustereddiffs = [PminTclustereddiffs metric_cluster_cifti_surfacearea(PminT,ID-.001,ID+.001,sizethresh,SAnames)];
%     end
    
    IminT = (InfomapData(:,s) - TemplateData(:,s)) .* (TemplateData(:,s) > 0) .* (InfomapData(:,s) > 0);
    IDs = unique(IminT); IDs(IDs==0) = [];
    IminTclustereddiffs = zeros(size(IminT,1),0);
    for ID = IDs(:)'
         IminTclustereddiffs = [IminTclustereddiffs metric_cluster_cifti_surfacearea(IminT,ID-.001,ID+.001,sizethresh,SAnames)];
    end
    
    MSCname = MSCnames{s};
    tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
    [subjects tmasks] = textread(tmaskfile,'%s %s');
    for i = 1:length(subjects)
        if i == 1
            data = cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
            tmask = load(tmasks{i});
        else
            data = [data cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'])];
            tmask = [tmask; load(tmasks{i})];
        end
    end
    
    data = data(1:ncortverts,logical(tmask));
    data(ncortverts+1:66697,:) = 0;
    
%     PminItcs = zeros(size(PminIclustereddiffs,2),size(data,2));
%     for i = 1:size(PminIclustereddiffs,2)
%         PminItcs(i,:) = mean(data(logical(PminIclustereddiffs(:,i)),:),1);
%         Pval = mode(PerceptronData(logical(PminIclustereddiffs(:,i)),s));
%         Ival = mode(InfomapData(logical(PminIclustereddiffs(:,i)),s));
%         PminIlabels{end+1} = ['Sub' num2str(s) '_P=' names{Pval} '_I=' names{Ival}];
%     end
%     map = FisherTransform(paircorr_mod(data',PminItcs'));
%     map(logical(PminIclustereddiffs)) = 0;
%     PminIcorrelmaps = [PminIcorrelmaps map];
%     
%     PminTtcs = zeros(size(PminTclustereddiffs,2),size(data,2));
%     for i = 1:size(PminTclustereddiffs,2)
%         PminTtcs(i,:) = mean(data(logical(PminTclustereddiffs(:,i)),:),1);
%         Pval = mode(PerceptronData(logical(PminTclustereddiffs(:,i)),s));
%         Tval = mode(TemplateData(logical(PminTclustereddiffs(:,i)),s));
%         PminTlabels{end+1} = ['Sub' num2str(s) '_P=' names{Pval} '_T=' names{Tval}];
%     end
%     map = FisherTransform(paircorr_mod(data',PminTtcs'));
%     map(logical(PminTclustereddiffs)) = 0;
%     PminTcorrelmaps = [PminTcorrelmaps map];
    
    IminTtcs = zeros(size(IminTclustereddiffs,2),size(data,2));
    for i = 1:size(IminTclustereddiffs,2)
        IminTtcs(i,:) = mean(data(logical(IminTclustereddiffs(:,i)),:),1);
        Ival = mode(InfomapData(logical(IminTclustereddiffs(:,i)),s));
        Tval = mode(TemplateData(logical(IminTclustereddiffs(:,i)),s));
        IminTlabels{end+1} = ['Sub' num2str(s) '_I=' names{Ival} '_T=' names{Tval}];
    end
    map = FisherTransform(paircorr_mod(data',IminTtcs'));
    map(logical(IminTclustereddiffs)) = 0;
    IminTcorrelmaps = [IminTcorrelmaps map];
end


%     cifti_write_wHDR(PminIcorrelmaps,[],['MSC_PminIcorrelmaps'])
%     fid = fopen('PminIlabels.txt','at');
%     fclose(fid);
%     for i = 1:length(PminIlabels)
%         dlmwrite('PminIlabels.txt',PminIlabels{i},'delimiter','','-append')
%     end
%     system('wb_command -cifti-convert-to-scalar MSC_PminIcorrelmaps.dtseries.nii ROW MSC_PminIcorrelmaps.dscalar.nii -name-file PminIlabels.txt');
%     
%     cifti_write_wHDR(PminTcorrelmaps,[],['MSC_PminTcorrelmaps'])
%     fid = fopen('PminTlabels.txt','at');
%     fclose(fid);
%     for i = 1:length(PminTlabels)
%         dlmwrite('PminTlabels.txt',PminTlabels{i},'delimiter','','-append')
%     end
%     system('wb_command -cifti-convert-to-scalar MSC_PminTcorrelmaps.dtseries.nii ROW MSC_PminTcorrelmaps.dscalar.nii -name-file PminTlabels.txt');
    
    cifti_write_wHDR(IminTcorrelmaps,[],['MSC_IminTcorrelmaps'])
    fid = fopen('IminTlabels.txt','at');
    fclose(fid);
    for i = 1:length(IminTlabels)
        dlmwrite('IminTlabels.txt',IminTlabels{i},'delimiter','','-append')
    end
    system('wb_command -cifti-convert-to-scalar MSC_IminTcorrelmaps.dtseries.nii ROW MSC_IminTcorrelmaps.dscalar.nii -name-file IminTlabels.txt');

    
    
    
    