
outfolder = ['/home/data/Analysis/MAV_analysis/convergence/'];
cd(outfolder);


initsubjects = textread('/home/data/Analysis/MAV_analysis/MAV_list.txt','%s');
initsubjects(strcmp(initsubjects,'MAV019')) = [];
initsubjects(strcmp(initsubjects,'MAV066')) = [];
initsubjects(strcmp(initsubjects,'MAV065')) = [];

[~,~,data]=xlsread('/home/data/Analysis/MAV_analysis/MAV_AssessmentData_cleanedv2.xlsx');
subjects = cell(0,1);
subs_withdata_indvec = false(size(data,1),1);
for s = 1:size(data,1)
    subs_withdata_indvec(s) = any(strcmp(data{s,1},initsubjects));
    if subs_withdata_indvec(s)
        subjects(end+1,1) = {data{s,1}};
    end
end


sessions_totest = [1 2 3];

tracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz');
%tracts = load_untouch_nii_2D('/home/data/Analysis/MAV_analysis/JHU-ICBM-labels-2mm-ero.nii.gz');
FA_similarity = ones(180,length(sessions_totest),length(subjects)) .* NaN;

for subnum = 1:length(subjects)
    
    subject = subjects{subnum};
    
    prevstring = [];
    
    WMmask = load_untouch_nii_2D(['/home/data/subjects/' subject '/DTI/aparc+aseg_cerebralwm_mask_222.nii.gz']);
    sub_tracts = tracts.img .* uint8(WMmask.img);
    tractIDs = unique(sub_tracts); tractIDs(tractIDs<1) = [];

    [runnames,~] = textread(['/home/data/subjects/' subject '/DTI/DTIruns_sessions.txt'],'%s%s');
    
    num_set_aside = floor(length(runnames) / 2);
    set_aside_combinations = nchoosek([1:length(runnames)],num_set_aside);

    
    FA_intracts_byrun = zeros(length(tractIDs),length(runnames));
    for r = 1:length(runnames)
        FA = load_untouch_nii_2D(['/home/data/subjects/' subject '/DTI/' runnames{r}]);
        for t = 1:length(tractIDs)
            FA_intracts_byrun(t,r) = mean(FA.img(sub_tracts==tractIDs(t)));
        end
    end
    
    for numsessionsind = 1:length(sessions_totest)
        if length(runnames) >= (num_set_aside + sessions_totest(numsessionsind))
        counter = 0;
        for combinationnum = 1:size(set_aside_combinations,1)
            set_aside_FA = mean(FA_intracts_byrun(:,set_aside_combinations(combinationnum,:)),2);
            
            all_test_samples = setdiff([1:length(runnames)],set_aside_combinations(combinationnum,:));
            test_samples = nchoosek(all_test_samples,sessions_totest(numsessionsind));
            for samplenum = 1:size(test_samples,1)
                counter = counter + 1;
                string = [subject ': sample ' num2str(counter) ', ' num2str(sessions_totest(numsessionsind)) ' sessions'];
                fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
                prevstring = string;
                test_FA = mean(FA_intracts_byrun(:,test_samples(samplenum,:)),2);
                FA_similarity(counter,numsessionsind,subnum) = paircorr_mod(test_FA,set_aside_FA);
            end
        end
        end
        save([outfolder '/FA_similarity_metrics.mat'],'FA_similarity')
    end
    disp(' ')
end
            
%%
outfolder = ['/home/data/Analysis/MAV_analysis/convergence/'];
FA_similarity = smartload([outfolder '/FA_similarity_metrics.mat']);


sessions_totest = [1 2 3];
outfolder = pwd;
colors = distinguishable_colors(length(subjects));

h = figure('Color','white','position',[100 100 2000 1200],'DefaultAxesFontSize',40);
hold on

final_similarities_toplot = zeros(length(subjects),length(sessions_totest));

for subnum = 1:length(subjects)
    legendnames{subnum} = subjects{subnum};
    meancorrmat = nanmean(FA_similarity(:,:,subnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(FA_similarity(:,:,subnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(FA_similarity,0,1);
    final_similarities_toplot(subnum,:) = meancorrmat;
    plot(sessions_totest,final_similarities_toplot(subnum,:),'-x','MarkerSize',20,'Color',colors(subnum,:),'LineWidth',5)
end
xlim([.9 3])
set(gca,'XTick',sessions_totest)
set(gca,'YTick',[.92 : .02 :1])
set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3,'FontName','Helvetica')
title(gca,['FA Reliability'])
xlabel('DTI Sessions','Fontweight','bold','FontSize',50,'FontName','Helvetica')
ylabel('Correlation to split half','Fontweight','bold','FontSize',50,'FontName','Helvetica')
%ylim([.4 .95])
%legend(legendnames,'Location','SouthEast')

export_fig(gca,[outfolder '/FA Similarity.pdf'])

