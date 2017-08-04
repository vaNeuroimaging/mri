
outfolder = ['/home/data/Analysis/MAV_analysis/convergence/'];

iterations = 1000;
%splithalf_quant = 70;
datalength_totest = [2.5 5 10 : 10 : 100];
TR = 3.0;



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

%splithalf_quant_frames = round(splithalf_quant * 60 / TR);
datalength_totest_frames = round(datalength_totest * 60 / TR);
%%
cd(outfolder);


corrmat_similarity = zeros(iterations,length(datalength_totest),length(subjects));

parcels_LR = '/home/data/atlases/Group_parcellation/Parcels_LR.dtseries.nii';
parcels_struct = ft_read_cifti_mod(parcels_LR);
parcels = parcels_struct.data;
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];
indmat = triu(true(length(parcelIDs)),1);

for subnum = 1:length(subjects)
    
    subject = subjects{subnum};
    
    tmask = logical(load(['/home/data/subjects/' subject '/fc_processed/RSFC_all_tmask.txt']));
    [runs,sessions] = textread(['/home/data/subjects/' subject '/fc_processed/RSFC_runs_sessions.txt'],'%f%f');
    runs = runs(tmask);
    data = ft_read_cifti_mod(['/home/data/subjects/' subject '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
    alltcs = zeros(size(data.data,2),length(parcelIDs));
    for IDnum = 1:length(parcelIDs)
        alltcs(:,IDnum) = mean(data.data(parcels==parcelIDs(IDnum),:),1);
    end
    clear data
    
    runnums = unique(runs);
    numsessions = length(runnums);
    numframes = zeros(numsessions,1);
%    tcs = cell(numsessions,1);
    for i = 1:numsessions
        numframes(i) = nnz(runs==runnums(i));
%        tcs{s} = alltcs(runs==runnums(i),:);
    end
    
    
    
    prevstring = [];
    
    for amountnum = 1:length(datalength_totest_frames)
        for iter = 1:iterations
            
            string = [subject ': iteration ' num2str(iter) ', ' num2str(datalength_totest(amountnum)) ' minutes of data'];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            
            randorder = randperm(numsessions);
            
            half_num_sessions = ceil(length(randorder)*1/2);
            
            set_aside_sessions = randorder(1:half_num_sessions);
            set_aside_runnums = runnums(set_aside_sessions);
            
            set_aside_inds = false(length(runs),1);
            for r = set_aside_runnums(:)'
                set_aside_inds = set_aside_inds | (runs==r);
            end
            
            
            
            all_test_sessions = setdiff(randorder,set_aside_sessions);
            all_test_sessions = all_test_sessions(randperm(length(all_test_sessions)));
            all_test_runnums = runnums(all_test_sessions);
            all_test_numframes = numframes(all_test_sessions);
            
            if sum(all_test_numframes) < datalength_totest_frames(amountnum) %nnz(all_test_inds)
                corrmat_similarity(iter,amountnum,subnum) = NaN;
            else
                test_inds = false(length(runs),1);
                framecount = 0;
                for r = all_test_runnums(:)'
                    if (framecount + numframes(runnums==r)) <= datalength_totest_frames(amountnum)
                        test_inds = test_inds | (runs==r);
                        framecount = framecount + numframes(runnums==r);
                    else
                        framesneeded = datalength_totest_frames(amountnum) - framecount;
                        cumulative_frames = cumsum(single(runs==r));
                        test_inds = test_inds | ((runs==r) & (cumulative_frames <= framesneeded));
                        break
                    end
                end
                
                
                
                set_aside_corrmat = paircorr_mod(alltcs(set_aside_inds,:));
                set_aside_corrmat(isnan(set_aside_corrmat)) = 0;
                set_aside_corrmat = FisherTransform(set_aside_corrmat);
                
                test_corrmat = paircorr_mod(alltcs(test_inds,:));
                test_corrmat(isnan(test_corrmat)) = 0;
                test_corrmat = FisherTransform(test_corrmat);
                
                corrmat_similarity(iter,amountnum,subnum) = paircorr_mod(set_aside_corrmat(indmat),test_corrmat(indmat));
                clear set_aside_corrmat test_corrmat
                
            end
            
            
%             iteration_failed = true;
%             
%             while iteration_failed
%                 
%                 randorder = randperm(numsessions);
%                 
%                 half_num_sessions = ceil(length(randorder)*1/2);
%                 
%                 %get set-aside data
%                 set_aside_sessions = randorder(1:half_num_sessions);
%                 set_aside_session_lengths = numframes(set_aside_sessions);
%                 [set_aside_session_lengths,sortorder] = sort(set_aside_session_lengths,'ascend');
%                 set_aside_sessions = set_aside_sessions(sortorder);
%                 
%                 if sum(set_aside_session_lengths) >= splithalf_quant_frames;
%                     iteration_failed = false;
%                     
%                     set_aside_tcs = zeros(0,length(parcelIDs));
%                     
%                     for sesscount = 1:length(set_aside_sessions);
%                         sessnum = set_aside_sessions(sesscount);
%                         left_tograb = splithalf_quant_frames - size(set_aside_tcs,1);
%                         thissess_length = size(tcs{sessnum},1);
%                         amount_tograb_fromsess = min([thissess_length (ceil(left_tograb ./ (length(set_aside_sessions) - sesscount + 1)))]);
%                         
%                         wiggle_room_insess = thissess_length - amount_tograb_fromsess;
%                         startpos = randi((wiggle_room_insess+1),1);
%                         set_aside_tcs((end+1):(end+amount_tograb_fromsess),:) = tcs{sessnum}(startpos : (amount_tograb_fromsess + startpos - 1),:);
%                     end
%                     
%                     
%                 end
%             end
%             
%             
%             %get test data
%             test_sessions = setdiff([1:length(subjectlist)],set_aside_sessions);
%             test_session_lengths = numframes(test_sessions);
%             [test_session_lengths,sortorder] = sort(test_session_lengths,'ascend');
%             test_sessions = test_sessions(sortorder);
%             
%             test_tcs = zeros(0,length(parcelIDs));
%             
%             for sesscount = 1:length(test_sessions);
%                 sessnum = test_sessions(sesscount);
%                 left_tograb = datalength_totest_frames(amountnum) - size(test_tcs,1);
%                 thissess_length = size(tcs{sessnum},1);
%                 amount_tograb_fromsess = min([thissess_length (ceil(left_tograb ./ (length(test_sessions) - sesscount + 1)))]);
%                 
%                 wiggle_room_insess = thissess_length - amount_tograb_fromsess;
%                 startpos = randi((wiggle_room_insess+1),1);
%                 test_tcs((end+1):(end+amount_tograb_fromsess),:) = tcs{sessnum}(startpos : (amount_tograb_fromsess + startpos - 1),:);
%             end
%             
%             if size(test_tcs,1) < datalength_totest_frames(amountnum);
%                 iteration_failed = true;
%             end
%             
%             
%             if iteration_failed
%                 corrmat_similarity(iter,amountnum,MSCnum) = NaN;
%                 
%                 
%             else
%                 
%                 
%                 
%                 
%                 
%                 set_aside_corrmat = paircorr_mod(set_aside_tcs);
%                 set_aside_corrmat(isnan(set_aside_corrmat)) = 0;
%                 set_aside_corrmat = FisherTransform(set_aside_corrmat);
%                 
%                 
%                 test_corrmat = paircorr_mod(test_tcs);
%                 test_corrmat(isnan(test_corrmat)) = 0;
%                 test_corrmat = FisherTransform(test_corrmat);
%                 
%                 
%                 
%                 corrmat_similarity(iter,amountnum,MSCnum) = paircorr_mod(set_aside_corrmat(indmat),test_corrmat(indmat));
%                 
%                 
%             end
            
            
            
            
        end
        
        save([outfolder '/similarity_metrics.mat'],'corrmat_similarity')
    end
    
    disp(' ')
    
end


%%
outfolder = ['/home/data/Analysis/MAV_analysis/convergence/'];
corrmat_similarity = smartload([outfolder '/similarity_metrics.mat']);


datalength_totest = [2.5 5 10 : 10 : 100];
outfolder = pwd;
colors = distinguishable_colors(length(subjects));

h = figure('Color','white','position',[100 100 2000 1200],'DefaultAxesFontSize',40);
hold on

final_similarities_toplot = zeros(length(subjects),length(datalength_totest));

for subnum = 1:length(subjects)
    legendnames{subnum} = subjects{subnum};
    meancorrmat = nanmean(corrmat_similarity(:,:,subnum),1);
    meancorrmat(meancorrmat==0) = NaN;
    meancorrmat(sum(isnan(corrmat_similarity(:,:,subnum)),1) > 900) = NaN;
    stdcorrmat = nanstd(corrmat_similarity,0,1);
    final_similarities_toplot(subnum,:) = meancorrmat;
    plot(datalength_totest,final_similarities_toplot(subnum,:),'-x','MarkerSize',20,'Color',colors(subnum,:),'LineWidth',5)
end

set(gca,'FontSize',40,'FontWeight','bold','LineWidth',3,'FontName','Helvetica')
title(gca,['RSFC Reliability'])
xlabel('Time (minutes)','Fontweight','bold','FontSize',50,'FontName','Helvetica')
ylabel('Correlation to split half','Fontweight','bold','FontSize',50,'FontName','Helvetica')
ylim([.4 1])
%legend(legendnames,'Location','SouthEast')

try export_fig(gca,[outfolder '/Correlation Similarity.pdf'])
catch
    savefig(gcf,[outfolder '/Correlation Similarity.fig'])
end

