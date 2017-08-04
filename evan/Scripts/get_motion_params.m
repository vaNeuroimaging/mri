subjects = '/home/data/subjects/processing_list_allcomplete.txt';
%subjects = '/home/data/subjects/list_forJDP.txt';
subjects = textread(subjects,'%s');

FDthresh = .2;
brainradius = 50;
TR = 2.2;
contiguoustimepoints = 5;

mean_FD = zeros(1,length(subjects));
pct_retained = zeros(1,length(subjects));
time_retained = zeros(1,length(subjects));

mean_FD_highdata = [];
pct_retained_highdata = [];
time_retained_highdata = [];

firstscans_totest = [1:10];
FD_testVall = cell(length(firstscans_totest),1);

for s = 1:length(subjects)
    subject = subjects{s};
    
%     preprocdir = ['/home/data/subjects/' subject '/preprocessed/'];
%     FD{s} = [];
%     tmask{s} = [];
%     
%     motionparamsfiles = dir([preprocdir '*.par']);
%     clear thisrun_FD
%     for filenum = 1:length(motionparamsfiles)
%         
%         thisrun_params = load([preprocdir motionparamsfiles(filenum).name]);
%         
%         thisrun_rot = thisrun_params(:,1:3);
%         thisrun_rot_mm = thisrun_rot * brainradius;
%         thisrun_params(:,1:3) = thisrun_rot_mm;
%         thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];
%         thisrun_FD{filenum} = [sum(abs(thisrun_params_delta),2)]';
%         
%         FD{s} = [FD{s} ; thisrun_FD{filenum}'];
%         
%         thisrun_tmask = (thisrun_FD{filenum}') < FDthresh;
%         
%         
%         chunkedges = [1; diff([thisrun_tmask;0])];
%         chunkstarts = find(chunkedges==1); chunkstops = (find(chunkedges==-1) -1);
%         chunksizes = zeros(length(chunkstarts),1);
%         chunkID = zeros(size(thisrun_tmask));
%         for chunknum = 1:length(chunkstarts)
%             chunkID(chunkstarts(chunknum) : chunkstops(chunknum)) = chunknum;
%             chunksizes(chunknum) = nnz(chunkID==chunknum);
%             if chunksizes(chunknum) < contiguoustimepoints
%                 thisrun_tmask(chunkID==chunknum) = 0;
%             end
%         end
%         
%         tmask{s} = [tmask{s} ; thisrun_tmask];
%     
%     end
%     
%     mean_FD(s) = mean(FD{s});
%     
%     
%     
%     pct_retained(s) = nnz(tmask{s}) / numel(tmask{s});
%     time_retained(s) = nnz(tmask{s}) * TR / 60;
% 
%     
%     if length(motionparamsfiles) > 40
%         mean_FD_highdata(end+1) = mean_FD(s);
%         pct_retained_highdata(end+1) = pct_retained(s);
%         time_retained_highdata(end+1) = time_retained(s);
%     end
%     
%     
%     for testnum = firstscans_totest
%         if (length(motionparamsfiles)/2) > testnum
%             test_FD = [thisrun_FD{1:testnum}];
%             test_pctretained = nnz(test_FD < FDthresh) / numel(test_FD);
%             FD_testVall{testnum}(end+1) = test_pctretained - pct_retained(s);
%         end
%         
%         test_names{testnum} = [num2str(testnum) ' session'];
%         
%     end
    
    


tmask{s} = load(['/home/data/subjects/' subject '/fc_processed/restingstate_all_tmask.txt']);
pct_retained(s) = nnz(tmask{s}) / numel(tmask{s});
time_retained(s) = nnz(tmask{s}) * TR / 60;

    
end

% figure; hist(mean_FD,10)
% figure; hist(pct_retained,10)
% figure; hist(time_retained,20)

% figure
% plotSpread(FD_testVall,'xNames',test_names,'distributionColors',repmat({'k'},1,length(test_names)),'MarkerSize',20)
% set(gca,'FontSize',20)
    
    