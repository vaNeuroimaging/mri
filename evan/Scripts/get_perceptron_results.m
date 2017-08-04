surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects, ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects, tmasks] = textread(tmaskfile,'%s %s');

subjects = subjects(81:100);%(1:120);%(121:end);
tmasks = tmasks(81:100);%(1:120);%(121:end);
ciftifiles = ciftifiles(81:100);%(1:120);%(121:end);

sizethreshmm = 100;

load 80train_20test_out.mat
load PCA_results_80.mat

consensus = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
IDs = unique(consensus); IDs(IDs==0) = [];

out = zeros(size(consensus,1),length(subjects));

for s = 1:length(subjects)
    disp(['Subject ' num2str(s)])
    tmask = load(tmasks{s});
    data = cifti_read(ciftifiles{s});
    data = data(:,logical(tmask))';
    
    corrmat = FisherTransform(paircorr_mod(data));
    corrmat(isnan(corrmat)) = 0;
    
    submapping = [ones(size(corrmat,1),1) corrmat * eigenvectors];
    submapping = submapping * weights1;
    submapping = tanh(.1 * submapping); submapping(:,1) = 1;
    submapping = submapping * out_weights1;
    
    for i = 1:size(submapping,2)
        [ncl(:,i),coutcl(:,i)]=hist(submapping(:,i),1000);
        r=range(coutcl(:,i));
        newcoutcl(:,i)=[min(coutcl(:,i))-.5*r;coutcl(:,i);max(coutcl(:,i))+.5*r];
        %     [ncl(:,i),coutclfix(:,i)]=hist(data(:,i),newcoutcl(:,i));
        %     data_cdf(:,i)=cumsum((ncl(:,i)))/sum((ncl(:,i)));
        %     data_U01(:,i)=interp1(coutclfix(:,i),data_cdf(:,i),data(:,i));
        [nclfix,coutclfix(:,i)]=hist(submapping(:,i),newcoutcl(:,i));
        data_cdf(:,i)=cumsum((nclfix))/sum((nclfix));
        submapping_uniform(:,i)=interp1(coutclfix(:,i),data_cdf(:,i),submapping(:,i));
    end
    clear ncl coutcl newcoutcl coutclfix data_cdf
    
    [ign maxi] = max(submapping_uniform,[],2);
    
    out(:,s) = IDs(maxi);
        
end

cifti_write_wHDR(out,[],'Perceptron_results_subs81thru100')
cleaned = Cifti_clean('Perceptron_results_subs81thru100.dtseries.nii',sizethreshmm);
cifti_write_wHDR(cleaned,[],['Perceptron_results_subs81thru100_clean' num2str(sizethreshmm) 'mm'])



