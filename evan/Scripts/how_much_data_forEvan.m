%% Setup Variables

[subjects tmasks] = textread('/data/cn3/steven/CoE/FCPROCESS_init/COHORTSELECT/NEW_TMASKLIST.txt','%s%s');
[subjects ciftifiles] = textread('/data/nil-bluearc/GMT/Evan/CoE/cifti_datalist.txt','%s%s');
outfolder = '/data/nil-bluearc/GMT/Evan/CoE/Analysis/Convergence/';
runlengths = 132; %frames
TR = 2.5;
parcelsfile = '/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii';
calc_timecourses = false;
iter = 100;

parcels = ft_read_cifti_mod(parcelsfile); parcels = parcels.data;
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];

runtime = runlengths * 2.5 / 60; %min

for subnum = 1:2%:length(subjects)
    
    %% Calculate timecourses
    disp([subjects{subnum} ': calculating parcel timecourses'])
    
    tmask = load(tmasks{subnum});
    nsess = length(tmask) / runlengths;
    
    if calc_timecourses
        subdata = ft_read_cifti_mod(ciftifiles{subnum}); subdata = subdata.data;
        timecourses = zeros(size(subdata,2),length(parcelIDs));
        
        for parcelnum = 1:length(parcelIDs)
            timecourses(:,parcelnum) = mean(subdata(parcels==parcelIDs(parcelnum),:),1);
            %tmask_bysess{sessnum} = tmask(sessinds);
        end
        
        
        for sessnum = 1:nsess
            sesind{sessnum} = [1:runlengths] + (runlengths * (sessnum-1));
        end
        
        clear subdata
        
        save([subjects{subnum} '_timecourses.mat'],'timecourses')
        save([subjects{subnum} '_sesind.mat'],'sesind')
        
    else
        timecourses = smartload([subjects{subnum} '_timecourses.mat']);
        sesind = smartload([subjects{subnum} '_sesind.mat']);
    end
    
    %
    %
    %
    %
    %
    %
    %
    %
    %
    % %% Load alreadly calculated timecourses
    % waterdir = '/data/cn4/laumannt/Poldrome/Poldrome_parcels';
    % [sessions tmasks] = textread('/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/poldrome_allses_selected_final_TMASKLIST.txt','%s%s');
    %
    % timecourses = smartload([waterdir '/parcel_timecourses_LR_mirpad_startpos50.mat']);
    
    
    %% Compare corrmat from subsets of sessions to corrmat from other half, concatenated, with replacement, random subjects for both halves, including same amount of time from more sesions
    %sessnum = length(sessions);
    maskmat = ones(size(timecourses,2));
    maskmat = triu(maskmat,1);
    %lengthtime = [518 259 172 129 103 86 74];
    factor = [1 2 4 8];
    %lengthtime = [518 259 129 64];
    halfsessionlength = runtime * floor(nsess/2);
    lengthtime = floor(runlengths ./ factor);
    totaltimeuse = floor(nsess ./ 2 ./ factor);
    
    
    clear corr_subset_to_cohhalf2
    %timecourses = [];
    %tmask = [];
    % for s = 1:nsess
    %     tmask = tmask_bysess{s};
    %     tmask_use = [tmask_use; tmask];
    %     time_concat = [time_concat; watertime_both{s}];
    %     sesind{s} = ((s-1)*518+1):s*518;
    % end
    count = 1;
    prevstring = [];
   
    
    for l = 1%:length(factor)
        
        for n = 1:iter
            string = [subjects{subnum} ': iterating data vs other half comparisons: iteration ' num2str(n)];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            
            cohhalf2 = randperm(nsess,floor(nsess/2))';
            indhalf2 = zeros(nsess,1);
            indhalf2(cohhalf2) = 1;
            indhalf1 = ~indhalf2;
            cohhalf1 = find(indhalf1);
            
            randsub = randperm(length(cohhalf1)); %Randomize order of sessions
            cohhalf1 = cohhalf1(randsub);
            lengthind = 1:lengthtime(1); % Full length for half comparing to
            % Calculate second half of data
            induse =[];
            for c = 1:length(cohhalf2)
                startpos = 0;
                subind = cohhalf2(c);
                induse = [induse sesind{subind}(lengthind+startpos)];
                %induse = [induse sesind{subind}];
                
            end
            time_use = timecourses(induse,:); % Use selected indices
            tmask_use_ind = tmask(induse); % Use selected indices
            %concat_corrmat_cohhalf2_startpos50(:,:,n) = FisherTransform(corrcoef(time_use(logical(tmask_use_ind),:)));
            concat_corrmat_cohhalf2 = FisherTransform(corrcoef(time_use(logical(tmask_use_ind),:)));
            %concat_corrmat_cohhalf2 = FisherTransform(corrcoef(time_concat(logical(tmask_use),:)));
            %cohhalf1 = randperm(subnum,subnum/2)';
            % Calculate correlation matrix sampled from first half of data
            for t = 1:totaltimeuse(l) %floor((subnum/2)/l)
                %disp(['time ' num2str(t*10)])
                lengthind = 1:lengthtime(l);
                %lengthind = 11:508; % Randomly sample
                %p = randperm(length(lengthind));
                %lengthind = lengthind(p(1:lengthtime(l)));
                %lengthind = (508-lengthtime(l)+1):508; % Use indices from back half of sessions
                induse = [];
                numsubuse = t*factor(l);
                for b = 1:numsubuse
                    subind = cohhalf1(b);
                    if l == 1
                        startpos = 0;
                    else
                        startpos = 0;
                        %   startpos = randperm(518-lengthtime(l)+1,1)-1; % Start from a random position in each session
                    end
                    induse = [induse sesind{subind}(lengthind+startpos)];
                    %induse = [induse sesind{subind}];
                    %sesind_step{b} = ((b-1)*length(stepind)+1):b*length(stepind);
                end
                time_use = timecourses(induse,:); % Use selected indices
                tmask_use_ind = tmask(induse); % Use selected indices
                concat_corrmat_sub_cohhalf1 = FisherTransform(corrcoef(time_use(logical(tmask_use_ind),:)));
                corr_subset_to_cohhalf2{count}(t,n) = paircorr_mod(concat_corrmat_sub_cohhalf1(logical(maskmat)),concat_corrmat_cohhalf2(logical(maskmat)));
                % concat_corrmat_sub_cohhalf1s{count}(:,:,t,n) = concat_corrmat_sub_cohhalf1;
            end
            
        end
        count = count + 1;
    end
    disp(' ')
    
    %% Plot correlation to other half with same amount of data with more vs. fewer sessions
    clear mean_cohhalf1_part_to_cohhalf2 mean_cohhalf1_full_to_cohhalf2
    clear std_cohhalf1_part_to_cohhalf2 std_cohhalf1_full_to_cohhalf2
    factor = [1 2 4 8];
    %lengthtime = [518 259 129 64];
    %lenghtime = [468 234 117 58];
    %totaltimeuse = [42 21 10 5];
    for l = 1%:length(lengthtime)
        mean_cohhalf1_part_to_cohhalf2{l} = mean(corr_subset_to_cohhalf2{l},2);
        std_cohhalf1_part_to_cohhalf2{l} = std(corr_subset_to_cohhalf2{l},[],2);
    end
    
    if subnum==1
    h = figure('Color','white','position',[1982 478 1352 504]);
    
    hold
    end
    linecolors = {'k','b','g','r','g','g','r'};
    time = runtime:runtime:floor(nsess/2)*runtime;
    plot(time,mean_cohhalf1_part_to_cohhalf2{1}(1:length(time)),'Color',linecolors{subnum},'LineWidth',3)
    %plot(time,mean_cohhalf1_part_to_cohhalf2{1}+std_cohhalf1_part_to_cohhalf2{1},'Color',linecolors{subnum},'LineStyle','--')
    %plot(time,mean_cohhalf1_part_to_cohhalf2{1}-std_cohhalf1_part_to_cohhalf2{1},'Color',linecolors{subnum},'LineStyle','--')
    %plot(time,mean_cohhalf1_part_to_cohhalf2{1}(1:length(time)),'Color','k','LineWidth',3)
    
    
    %  for l = [2 3 4] %[2 4 6 7]% 51:length(lengthtime)
    %
    %      time = runtime:runtime:totaltimeuse(l)*runtime;%floor((length(subjects)/2)/(l))*10;
    %      plot(time,mean_cohhalf1_part_to_cohhalf2{l},'--','Color',linecolors{l},'LineWidth',4)
    % %     %shadedErrorBar(time,cohhalf1_part_to_cohhalf2{l},{@mean,@std},linecolors{l})
    % %   %  plot(time,mean_cohhalf1_part_to_cohhalf2{l}+2*std_cohhalf1_part_to_cohhalf2{l},'Color',linecolors{l},'LineStyle','--')
    % %     %fill([time time],[(mean_cohhalf1_part_to_cohhalf2{l}+2*std_cohhalf1_part_to_cohhalf2{l}) (mean_cohhalf1_part_to_cohhalf2{l}-2*std_cohhalf1_part_to_cohhalf2{l})],linecolors{l})
    % %   % plot(time,mean_cohhalf1_part_to_cohhalf2{l}-2*std_cohhalf1_part_to_cohhalf2{l},'Color',linecolors{l},'LineStyle','--')
    % %
    %  end
    %  hleg = legend('10 min','5 min','2.5 min','1.25 min','location','southeast') %1.66 min','1.43 min')
    
    ylim([.5 1])
    xlim([0 150])
    
    set(gca,'FontSize',14,'FontWeight','bold','LineWidth',2)
    xlabel('Time (minutes)','Fontweight','bold','FontSize',14)
    ylabel('Corr to other half','Fontweight','bold','FontSize',14)
    
    %export_fig(gca,[outfolder '/' subjects{subnum} '_time_to_convergence.pdf'])
    
end
