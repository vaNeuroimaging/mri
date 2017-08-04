
outputfilenames = {'DA_lOT_eroded_connectivity.func.gii','FPC_lOT_eroded_connectivity.func.gii'};

ROIfiles = {'120consensus_DA_clusters_5_to_5_minsize20_eroded2.func.gii','120consensus_FPC_clusters_6_to_6_minsize20_eroded2.func.gii'};
ROIindices = [3 6];

for roinum = 1:length(ROIfiles)

    ROIs = gifti(ROIfiles{roinum});
    ROI{roinum} = ROIs.cdata(:,ROIindices(roinum));
end

HEMS = {'L';'R'};
hem = 1;

surftimedir = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort/';

cohortdir = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort/'];
[subjects tmasks] = textread([cohortdir '/NEW_nokids_TMASKLIST.txt'],'%q %q');

for s = 1:length(subjects)
    
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        tic
        timename = [subjects{s} '_BOLD_' HEMS{hem} '_smooth2.55_32k_fsLR'];
%         system(['caret_command64 -file-convert -format-convert ASCII ' timename '.func.gii'])
%         system(['awk ''NF > 25'' ' timename '.func.gii > ' timename '_noHEAD.func.gii'])
        tmask = load(tmasks{s});
        surf_timecourses = load([surftimedir timename '_noHEAD.func.gii']);
        surf_timecourses(:,1) = [];
        surf_timecourses = surf_timecourses(:,logical(tmask));
        
        for roinum = 1:length(ROIfiles)
        
            roi_timecourse = mean(surf_timecourses(logical(ROI{roinum}),:),1);
        
            roi_pattern{roinum}(:,s) = paircorr_mod(roi_timecourse', surf_timecourses');
        end
        
end

for roinum = 1:length(ROIfiles)

Fisher_roi_pattern = FisherTransform(roi_pattern{roinum});

mean_roi_pattern = mean(Fisher_roi_pattern,2);

save(gifti(single(mean_roi_pattern)),outputfilenames{roinum});
end