%%
trainingnums = [1:80];
testnums = [81:100];

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects, ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects, tmasks] = textread(tmaskfile,'%s %s');

subjects = subjects(trainingnums);%(81:100);%(121:end);
tmasks = tmasks(trainingnums);%(81:100);%(121:end);
ciftifiles = ciftifiles(trainingnums);%(81:100);%(121:end);


% consensus = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
% IDs = unique(consensus); IDs(IDs==0) = [];
% 
% ROIs = zeros(size(consensus,1),0);
% 
% ROI_IDs = zeros(0,length(IDs));
% 
% for IDnum = 1:length(IDs)
%     ID = IDs(IDnum);
%     IDclusters = metric_cluster_cifti(consensus,ID-.01,ID+.01,5);
%     ROIs(:,end+1 : end+size(IDclusters,2)) = IDclusters;
%     ROI_IDs(end+1 : end+size(IDclusters,2) , IDnum) = 1;
% end


Parcels = cifti_read('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii'); Parcels(end+1:66697) = 0;
Communities = cifti_read('/data/cn4/evan/Published_parcels/Parcel_Communities.dtseries.nii');
IDs = unique(Communities); IDs(IDs<=0) = [];

ROIs = zeros(size(Parcels,1),0);
%ROI_IDs = [];

for IDnum = 1:length(IDs)
    ID = IDs(IDnum);
    parcelnums_thisID = unique(Parcels(Communities==ID));
    
    for parcelnum = parcelnums_thisID(:)'
        ROIs(:,end+1) = single(logical(Parcels==parcelnum));
        %ROI_IDs(end+1,1) = IDnum;
    end
end



ROISubVertMat = zeros(size(ROIs,2) * length(subjects) , size(Parcels,1));

for s = 1:length(subjects)
    disp(['Subject ' num2str(s)])
    tmask = load(tmasks{s});
    data = cifti_read(ciftifiles{s});
    data = data(:,logical(tmask))';
    
    for ROInum = 1:size(ROIs,2)
        matrixindex = ((s-1) * size(ROIs,2)) + ROInum;
        ROItimecourse = mean(data(:,logical(ROIs(:,ROInum))),2);
        ROIcorr = paircorr_mod(ROItimecourse , data);
        ROIcorr(isnan(ROIcorr)) = 0;
        ROISubVertMat(matrixindex,:) = FisherTransform(ROIcorr);
    end
end

%save('Test20_ROISubVertMat.mat','ROISubVertMat','-v7.3')

%ROISub_IDmat = repmat(ROI_IDs,length(subjects),1);

%save('ROISub_IDmat.mat','ROISub_IDmat');
save(['ROISubVertMat_' num2str(length(trainingnums)) '.mat'],'ROISubVertMat','-v7.3')

[components, eigenvectors, perc_variance_explained] = PCA_reduction_specify_numcomps(ROISubVertMat,2500);

save(['PCA_results_' num2str(length(trainingnums)) '.mat'],'components', 'eigenvectors', 'perc_variance_explained','-v7.3')



%%

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects, ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects, tmasks] = textread(tmaskfile,'%s %s');

subjects = subjects(testnums);
tmasks = tmasks(testnums);
ciftifiles = ciftifiles(testnums);

Parcels = cifti_read('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii'); Parcels(end+1:66697) = 0;
Communities = cifti_read('/data/cn4/evan/Published_parcels/Parcel_Communities.dtseries.nii');
IDs = unique(Communities); IDs(IDs<=0) = [];

ROIs = zeros(size(Parcels,1),0);
ROI_IDs = [];

for IDnum = 1:length(IDs)
    ID = IDs(IDnum);
    parcelnums_thisID = unique(Parcels(Communities==ID));
    
    for parcelnum = parcelnums_thisID(:)'
        ROIs(:,end+1) = single(logical(Parcels==parcelnum));
        ROI_IDs(end+1,1) = IDnum;
    end
end

ROISubVertMat = zeros(size(ROIs,2) * length(subjects) , size(Parcels,1));

for s = 1:length(subjects)
    disp(['Subject ' num2str(s)])
    tmask = load(tmasks{s});
    data = cifti_read(ciftifiles{s});
    data = data(:,logical(tmask))';
    
    for ROInum = 1:size(ROIs,2)
        matrixindex = ((s-1) * size(ROIs,2)) + ROInum;
        ROItimecourse = mean(data(:,logical(ROIs(:,ROInum))),2);
        ROIcorr = paircorr_mod(ROItimecourse , data);
        ROIcorr(isnan(ROIcorr)) = 0;
        ROISubVertMat(matrixindex,:) = FisherTransform(ROIcorr);
    end
end

load(['PCA_results_' num2str(length(trainingnums)) '.mat'])

training_in = components';
test_in = ROISubVertMat * eigenvectors;
desired_train = repmat(ROI_IDs,length(trainingnums),1);
desired_test = repmat(ROI_IDs,length(testnums),1);

save([num2str(length(trainingnums)) 'train_' num2str(length(testnums)) 'test_in.mat'],'training_in','test_in','desired_test','desired_train')


MLP_Training([num2str(length(trainingnums)) 'train_' num2str(length(testnums)) 'test_in.mat'])
