surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects, ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects, tmasks] = textread(tmaskfile,'%s %s');

subjects = subjects(81:100);%(1:120);%(121:end);
tmasks = tmasks(81:100);%(1:120);%(121:end);
ciftifiles = ciftifiles(81:100);%(1:120);%(121:end);

sizethresh = 300;

PerceptronData = cifti_read('/data/cn4/evan/RestingState/Ind_variability/Perceptron/Perceptron_results_subs81thru100_clean100mm.dtseries.nii');
InfomapData = cifti_read('/data/cn4/evan/RestingState/Ind_variability/120_infomaps_0002_to_003_in_0002_xd30_BI/Consensus_Powercolors_cleaned_allsubs_rawassn_minsize400_regularize_cleaned.dtseries.nii');
InfomapData = InfomapData(:,81:100);
TemplateData = cifti_read('/data/cn4/evan/RestingState/Ind_variability/120_108_clean100mm/Templatematch_dice_bysubject.dtseries.nii');
TemplateData = TemplateData(:,81:100);

for s = 1%:length(subjects)
    disp(['Subject ' num2str(s)])
    PminI = (PerceptronData(:,s) - InfomapData(:,s)) .* (PerceptronData(:,s) > 0) .* (InfomapData(:,s) > 0);
    IDs = unique(PminI); IDs(IDs==0) = [];
    PminIclustereddiffs = zeros(size(PminI,1),0);
    for ID = IDs(:)'
         PminIclustereddiffs = [PminIclustereddiffs metric_cluster_cifti_surfacearea(PminI,ID-.001,ID+.001,sizethresh) ];
    end
    
    PminT = (PerceptronData(:,s) - TemplateData(:,s)) .* (PerceptronData(:,s) > 0) .* (TemplateData(:,s) > 0);
    IDs = unique(PminT); IDs(IDs==0) = [];
    PminTclustereddiffs = zeros(size(PminT,1),0);
    for ID = IDs(:)'
         PminTclustereddiffs = [PminTclustereddiffs metric_cluster_cifti_surfacearea(PminT,ID-.001,ID+.001,sizethresh)];
    end
    
    IminT = (InfomapData(:,s) - TemplateData(:,s)) .* (TemplateData(:,s) > 0) .* (InfomapData(:,s) > 0);
    IDs = unique(IminT); IDs(IDs==0) = [];
    IminTclustereddiffs = zeros(size(IminT,1),0);
    for ID = IDs(:)'
         IminTclustereddiffs = [IminTclustereddiffs metric_cluster_cifti_surfacearea(IminT,ID-.001,ID+.001,sizethresh)];
    end
    
    tmask = load(tmasks{s});
    data = cifti_read(ciftifiles{s});
    data = data(:,logical(tmask));
    
    PminItcs = zeros(size(PminIclustereddiffs,2),size(data,2));
    for i = 1:size(PminIclustereddiffs,2)
        PminItcs(i,:) = mean(data(logical(PminIclustereddiffs(:,i)),:),1);
    end
    PminIcorrelmaps = FisherTransform(paircorr_mod(data',PminItcs'));
    
    PminTtcs = zeros(size(PminTclustereddiffs,2),size(data,2));
    for i = 1:size(PminTclustereddiffs,2)
        PminTtcs(i,:) = mean(data(logical(PminTclustereddiffs(:,i)),:),1);
    end
    PminTcorrelmaps = FisherTransform(paircorr_mod(data',PminTtcs'));
    
    IminTtcs = zeros(size(IminTclustereddiffs,2),size(data,2));
    for i = 1:size(IminTclustereddiffs,2)
        IminTtcs(i,:) = mean(data(logical(IminTclustereddiffs(:,i)),:),1);
    end
    IminTcorrelmaps = FisherTransform(paircorr_mod(data',IminTtcs'));
    
    cifti_write_wHDR(PminIclustereddiffs,[],['Subject' num2str(s) '_PminIdiffs'])
    cifti_write_wHDR(PminIcorrelmaps,[],['Subject' num2str(s) '_PminIcorrelmaps'])
    
    cifti_write_wHDR(PminTclustereddiffs,[],['Subject' num2str(s) '_PminTdiffs'])
    cifti_write_wHDR(PminTcorrelmaps,[],['Subject' num2str(s) '_PminTcorrelmaps'])
    
    cifti_write_wHDR(IminTclustereddiffs,[],['Subject' num2str(s) '_IminTdiffs'])
    cifti_write_wHDR(IminTcorrelmaps,[],['Subject' num2str(s) '_IminTcorrelmaps'])
    
end
    
    
    
    