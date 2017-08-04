ncortverts = 59412;

MSCnames = {'MSC01','MSC02','MSC03','MSC04'};

sizethreshmm = 50;

load 120train_4test_kden_out.mat
load PCA_results_120_kden.mat

consensus = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
IDs = unique(consensus); IDs(IDs==0) = [];

out = zeros(size(consensus,1),length(MSCnames));

for s = 1:length(MSCnames)
    disp(['Subject ' num2str(s)])
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
    
    data = data(1:ncortverts,logical(tmask))';
    
    corrmat = FisherTransform(paircorr_mod(data));
    corrmat(isnan(corrmat)) = 0;
    
    submapping = [ones(size(corrmat,1),1) corrmat * eigenvectors];
    submapping = submapping * weights1;
    submapping = tanh(.1 * submapping); submapping(:,1) = 1;
    submapping = submapping * out_weights1;
    
%     for i = 1:size(submapping,2)
%         [ncl(:,i),coutcl(:,i)]=hist(submapping(:,i),1000);
%         r=range(coutcl(:,i));
%         newcoutcl(:,i)=[min(coutcl(:,i))-.5*r;coutcl(:,i);max(coutcl(:,i))+.5*r];
%         %     [ncl(:,i),coutclfix(:,i)]=hist(data(:,i),newcoutcl(:,i));
%         %     data_cdf(:,i)=cumsum((ncl(:,i)))/sum((ncl(:,i)));
%         %     data_U01(:,i)=interp1(coutclfix(:,i),data_cdf(:,i),data(:,i));
%         [nclfix,coutclfix(:,i)]=hist(submapping(:,i),newcoutcl(:,i));
%         data_cdf(:,i)=cumsum((nclfix))/sum((nclfix));
%         submapping_uniform(:,i)=interp1(coutclfix(:,i),data_cdf(:,i),submapping(:,i));
%     end
%     clear ncl coutcl newcoutcl coutclfix data_cdf
%     
%     [ign maxi] = max(submapping_uniform,[],2);
    submapping_transformed = 1./(1 + exp(-submapping));
    [ign maxi] = max(submapping_transformed,[],2);
    out(1:ncortverts,s) = IDs(maxi);
        
end

cifti_write_wHDR(out,[],'Perceptron_kden_results_MSC')
cleaned = Cifti_clean('Perceptron_kden_results_MSC.dtseries.nii',sizethreshmm);
cifti_write_wHDR(cleaned,[],['Perceptron_kden_results_MSC_clean' num2str(sizethreshmm) 'mm'])



