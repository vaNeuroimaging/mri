%% Calculate ANOVA of correlation maps over sessions, single hemisphere
% TR = 2.5;
% HEMS = {'L';'R'};
% hemname = {'LEFT';'RIGHT'};
% hemnamelow = {'left';'right'};

hem = 'L';
mintimepoints = 245;
numparts = 7;



wall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']);
wall = wall.cdata;
corticalindices = find(wall==0);

smooth = 2.55;

outputdir = '/data/cn4/evan/RestingState/Ind_variability/Ftest/';

cohortdir = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort'];
[subjects tmasks] = textread([cohortdir '/NEW_nokids_TMASKLIST.txt'],'%q %q');

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
ciftidir = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/';

SSE = zeros(1,length(corticalindices));
SSA = zeros(1,length(corticalindices));

%if strcmp(hem,'L')
    ciftiindices = 1:length(corticalindices);
%elseif strcmp(hem,'R')
%    ciftiindices = [29697:29696+length(corticalindices)];
%end


%parts = {1:259;260:518};
%parts = {1:74;75:148;149:222;223:296;297:370;371:444;445:518};

avg_corr_global = zeros(length(corticalindices));
goodsubs = [];
for s = 1:length(subjects)
    %tic
    disp(['Processing subject #' num2str(s)])
    tmask = load(tmasks{s});
    
    numtimepoints = nnz(tmask);
    
    if numtimepoints < mintimepoints
        
        disp('Not enough data. Skipping.')
        
    else
    
    goodsubs(end+1) = s;
        
    tp_per_part = floor(numtimepoints/numparts);
    
    timename = [ciftidir '/' subjects{s} '_BOLD_' hem '_surf_subcort_contracbll_smooth2.55_32k_fsLR'];
    evalc(['!' workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' timename '.dtseries.nii /data/cn4/evan/Temp/Temp.func.gii']);
    timecourse = gifti('/data/cn4/evan/Temp/Temp.func.gii');
    timecourse = timecourse.cdata';
    timecourse = timecourse(logical(tmask),ciftiindices);
    
    for p = 1:numparts
        
        parts{p,1} = [(tp_per_part*(p-1) + 1) : tp_per_part*p];
        timecourse_part = timecourse(parts{p},:);
        
        %time_temp = timecourse(parts{p},1:29696);
        %timecourse_part = time_temp(logical(tmask(parts{p})),:);
        time_corr_part = single(FisherTransform(corrcoef(timecourse_part)));
        time_corr_part(isnan(time_corr_part)) = 0;
        avg_corr_global = avg_corr_global + time_corr_part;
    end
    %system(['rm ' timename '.func.*']);
    %toc
    end
end

avg_corr_global = avg_corr_global./(numparts*length(goodsubs));
%%
for s_ind = 1:length(goodsubs)
    s = goodsubs(s_ind);
    avg_corr_winses = zeros(length(corticalindices));
    smooth = 2.55;
    %tic
    disp(['Processing subject #' num2str(s)])
    tmask = load(tmasks{s});
    
    numtimepoints = nnz(tmask);
    tp_per_part = floor(numtimepoints/numparts);
    
    timename = [ciftidir '/' subjects{s} '_BOLD_LR_surf_subcort_smooth' num2str(smooth) '_32k_fsLR'];
    system([workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' timename '.dtseries.nii /data/cn4/evan/Temp/Temp.func.gii'])
    timecourse = gifti('/data/cn4/evan/Temp/Temp.func.gii');
    timecourse = timecourse.cdata';
    timecourse = timecourse(logical(tmask),ciftiindices);
    
    for p = 1:numparts
        
        parts{p,1} = [(tp_per_part*(p-1) +1) : tp_per_part*p];
        timecourse_part = timecourse(parts{p},:);
        
        %time_temp = timecourse(parts{p},1:29696);
        %timecourse_part = time_temp(logical(tmask(parts{p})),:);
        time_corr_part = single(FisherTransform(corrcoef(timecourse_part)));
        time_corr_part(isnan(time_corr_part)) = 0;
        avg_corr_winses = avg_corr_winses + time_corr_part;
    end
    avg_corr_winses = avg_corr_winses./numparts;
    SSA = SSA + numparts*sum((avg_corr_winses-avg_corr_global).^2);
    for p = 1:numparts
        
        timecourse_part = timecourse(parts{p},:);
        
        %time_temp = timecourse(parts{p},1:29696);
        %timecourse_part = time_temp(logical(tmask(parts{p})),:);
        time_corr_part = single(FisherTransform(corrcoef(timecourse_part)));
        time_corr_part(isnan(time_corr_part)) = 0;
        SSE = SSE + sum((time_corr_part-avg_corr_winses).^2);
        
    end
    %system(['rm ' timename '.func.*']);
    %toc
end
clear avg_corr_winses
%save([outputdir '/SSE_7part_' hem '.mat'],'SSE','-v7.3')
%save([outputdir '/SSA_7part_' hem '.mat'],'SSA','-v7.3')


MSA = SSA./(length(goodsubs)-1);
MSE = SSE./(length(goodsubs).*(numparts-1));
F = MSA./MSE; 

%save([outputdir '/F_7part_L.mat'],'F','-v7.3')


metric = zeros(size(wall));
metric(corticalindices,1) = F;

save(gifti(single(metric)),[outputdir '/F_7part_' hem '.func.gii'])