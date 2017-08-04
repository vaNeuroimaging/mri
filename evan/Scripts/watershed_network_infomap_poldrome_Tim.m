%% Load edge, upper complete, minima detection, watershed
%dir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients_cifti/all48';
%outputdir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients_cifti/all48/watershed';
minimadist = 5;
HEM = {'L';'R'};
for c = 1:3
    dir = ['/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/C' num2str(c)];
    outputdir = ['/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/C' num2str(c) '/watershed'];
    
    for h = 1:2
        edge=gifti([dir '/avgcorrofcorr_smooth2.55_allgrad_' HEM{h} '_smooth2.55_edge_avg.func.gii']);
        %edge=gifti([dir '/avgcorrofcorr_allgrad_' HEM(h) '_smooth2.55_edges_avg.func.gii']);
        edge = edge.cdata;
        [fUC fUC_smooth] = upper_completion(edge);
        cd(dir)
        save(gifti(single(fUC_smooth)),[dir '/avgcorrofcorr_smooth255_allgrad_' HEM{h} '_smooth255_edge_avg_uc_smooth.func.gii']);
        %save(gifti(single(fUC_smooth)),[outputdir '/avgcorrofcorr_allgrad_' HEM(h) '_smooth2.55_edges_avg_uc_smooth.func.gii']);
        minimametric = metric_minima(fUC_smooth,minimadist);
        cd(dir)
        %save(gifti(single(minimametric)),[dir '/avgcorrofcorr_smooth255_allgrad_' HEM '_smooth255_edge_avg_uc_smooth_minima' num2str(minimadist) '.func.gii']);
        %label = watershed_algorithm(fUC_smooth,minimametric,200,1,dir,['avgcorrofcorr_smooth255_allgrad_' HEM '_smooth255_edge_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_']);
        save(gifti(single(minimametric)),[outputdir '/avgcorrofcorr_smooth255_allgrad_' HEM{h} '_smooth255_edges_avg_uc_smooth_minima' num2str(minimadist) '.func.gii']);
        label = watershed_algorithm(fUC_smooth,minimametric,200,1,outputdir,['avgcorrofcorr_smooth255_allgrad_' HEM{h} '_smooth255_edges_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_']);
    end
end
%% Watershed-based networks, Left and Right
cohs = {'C1';'C2';'C3'};
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
waterdir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/all48/watershed';
outputdir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/all48/watershed/network';

surftimedir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/all48';
minimadist = 4;

cd(surftimedir)
[subjects tmasks] = textread('/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS/COHORTSELECT/NEW_TMASKLIST.txt','%s%s');

clear water_corrmat
for s = 1:length(subjects)
    clear watertime watertime_both
    for hem = 1:2
        watershed = gifti([waterdir '/avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_watershed.func.gii']);
        watershed = watershed.cdata;
        maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
        mask = gifti(maskname);
        mask = mask.cdata;
        watershed(logical(mask)) = 0;
        waternum = unique(watershed);
        waternum(waternum==0) = [];
        
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        tic
        timename = [subjects{s} '_BOLD_' HEMS{hem} '_smooth2.55_32k_fsLR'];
        system(['caret_command64 -file-convert -format-convert ASCII ' timename '.func.gii'])
        system(['awk ''NF > 25'' ' timename '.func.gii > ' timename '_noHEAD.func.gii'])
        tmask = load(tmasks{s});
        surf_timecourses = load([timename '_noHEAD.func.gii']);
        surf_timecourses(:,1) = [];
        surf_timecourses = surf_timecourses(:,logical(tmask));
        
        for i = 1:length(waternum)
            
            waterind = find(watershed==waternum(i));
            watertime{hem}(:,i) = mean(surf_timecourses(waterind,:))';
        end
        
        
    end
    
    watertime_both = [watertime{1} watertime{2}];
    water_corrmat(:,:,s) = paircorr_mod(watertime_both);
    system(['rm ' timename '_noHEAD.func.gii'])
end
   


%all_water_corrmat = reshape(water_corrmat,[size(watertime_both,2) size(watertime_both,2) 120]);
all_water_corrmat = water_corrmat;

save([outputdir '/avgcorrofcorr_allgrad_LR_smooth255_edges_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_watershed_corrmat.mat'],'all_water_corrmat')

mean_all_water_corrmat = mean(all_water_corrmat,3);
save([outputdir '/avgcorrofcorr_allgrad_LR_smooth255_edges_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_watershed_meancorrmat.mat'],'mean_all_water_corrmat')


%% Watershed-based networks, Left and Right, corrmap-based
cohs = {'C1';'C2';'C3'};
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
waterdir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/all48/watershed';
outputdir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/all48/watershed/network';

surftimedir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/all48';

cd(surftimedir)
[subjects tmasks] = textread('/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS/COHORTSELECT/NEW_TMASKLIST.txt','%s%s');

clear water_corrmat
for s = 1:length(subjects)
    clear watertime watertime_both
    for hem = 1:2
        watershed = gifti([waterdir '/avgcorrofcorr_smooth2.55_allgrad_' HEMS{hem} '_smooth2.55_edge_avg_uc_smooth_minima5_iter200_frach1_watershed.func.gii']);
        watershed = watershed.cdata;
        maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
        mask = gifti(maskname);
        mask = mask.cdata;
        watershed(logical(mask)) = 0;
        waternum = unique(watershed);
        waternum(waternum==0) = [];
        
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        tic
        timename = [subjects{s} '_BOLD_' HEMS{hem} '_smooth2.55_32k_fsLR'];
        system(['caret_command64 -file-convert -format-convert ASCII ' timename '.func.gii'])
        system(['awk ''NF > 25'' ' timename '.func.gii > ' timename '_noHEAD.func.gii'])
        tmask = load(tmasks{s});
        surf_timecourses = load([timename '_noHEAD.func.gii']);
        surf_timecourses(:,1) = [];
        surf_timecourses = surf_timecourses(:,logical(tmask));
        
        for i = 1:length(waternum)
            
            waterind = find(watershed==waternum(i));
            watertime{hem}(:,i) = mean(surf_timecourses(waterind,:))';
        end
        
        
    end
    
    watertime_both = [watertime{1} watertime{2}];
    water_corrmat(:,:,s) = paircorr_mod(watertime_both);
    system(['rm ' timename '_noHEAD.func.gii'])
end
   


%all_water_corrmat = reshape(water_corrmat,[size(watertime_both,2) size(watertime_both,2) 120]);
all_water_corrmat = water_corrmat;

save([outputdir '/avgcorrofcorr_allgrad_LR_smooth2.55_edges_avg_uc_smooth_minima5_iter200_frach1_watershed_corrmat.mat'],'all_water_corrmat')

mean_all_water_corrmat = mean(all_water_corrmat,3);
save([outputdir '/avgcorrofcorr_allgrad_LR_smooth2.55_edges_avg_uc_smooth_minima5_iter200_frach1_watershed_meancorrmat.mat'],'mean_all_water_corrmat')
%% Minima-ROI-based networks
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
waterdir = '/data/cn4/laumannt/watershed_network';
% 
% maskname = '/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii';
% 
% mask = gifti(maskname);
% mask = mask.cdata;
% % 
waterdir = '/data/cn4/laumannt/watershed_network/FCPROCESS_NEW/';

surftimedir = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort';

cd(surftimedir)
    
cohortdir = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort'];
[subjects tmasks] = textread([cohortdir '/NEW_nokids_TMASKLIST.txt'],'%q %q');
clear minima_corrmat 
for s = 1:length(subjects)
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        clear meansurf_timecourse
        for hem = 1:2
            minima_rois = gifti([waterdir '/avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center_rois.func.gii']);
            minima_rois = minima_rois.cdata;
%         maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
%         mask = gifti(maskname);
%         mask = mask.cdata;
%         watershed(logical(mask)) = 0;
%         waternum = unique(watershed);
%         waternum(waternum==0) = [];
            timename = [subjects{s} '_BOLD_' HEMS{hem} '_smooth2.55_32k_fsLR'];
            system(['caret_command64 -file-convert -format-convert ASCII ' timename '.func.gii'])
            system(['awk ''NF > 25'' ' timename '.func.gii > ' timename '_noHEAD.func.gii'])
            tmask = load(tmasks{s});
            surf_timecourses = load([timename '_noHEAD.func.gii']);
            surf_timecourses(:,1) = [];
            surf_timecourses = surf_timecourses(:,logical(tmask));
            
            
            for n = 1:size(minima_rois,2)
                
                meansurf_timecourse{hem}(:,n) = mean(surf_timecourses(logical(minima_rois(:,n)),:));
            end
            
        end
            minima_timecourse_both = [meansurf_timecourse{1} meansurf_timecourse{2}];
            minima_corrmat(:,:,s) = paircorr_mod(minima_timecourse_both);
            
            system(['rm ' timename '_noHEAD.func.gii'])

        
end

all_minima_corrmat = minima_corrmat;

cd('/data/cn4/laumannt/watershed_network/FCPROCESS_NEW/')
save('avgcorrofcorr_smooth255_allgrad_LR_smooth255_edge_avg_uc_smooth_minima5_parcel_center_rois_corrmat.mat','all_minima_corrmat')

mean_all_minima_corrmat = mean(all_minima_corrmat,3);
save('avgcorrofcorr_smooth255_allgrad_LR_smooth255_edge_avg_uc_smooth_minima5_parcel_center_rois_meancorrmat.mat','mean_all_minima_corrmat')

%% Run Infomap
graphcluster('prmfile.txt','thr',[] ,1,0,'infomap');

%% Modify color assignments with miminum network size criterion
simple_assigns = modify_clrfile('simplify','rawassn.txt',4);

%% Assign communities to watersheds, LR together
waterdir = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED/gradients/all48/watershed';
minimadist = 4;
watershed = gifti([waterdir '/avgcorrofcorr_smooth255_allgrad_L_smooth255_edge_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_watershed.func.gii']);
watershed = watershed.cdata;
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
watershed(logical(mask)) = 0;
waternum_L = unique(watershed);
waternum_L(waternum_L==0) = [];

water_com_assigns = zeros(32492,size(simple_assigns,2));

for t = 1:size(simple_assigns,2)
    
    for i = 1:length(waternum_L)
        
        water_com_assigns(watershed==waternum_L(i),t) = simple_assigns(i,t);
    end
end
save(gifti(water_com_assigns),[waterdir '/avgcorrofcorr_allgrad_L_smooth255_edges_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_watershed_rawassn_minsize4.func.gii'])
                 

watershed = gifti([waterdir '/avgcorrofcorr_smooth255_allgrad_R_smooth255_edge_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_watershed.func.gii']);
watershed = watershed.cdata;
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
watershed(logical(mask)) = 0;
waternum_R = unique(watershed);
waternum_R(waternum_R==0) = [];

water_com_assigns = zeros(32492,size(simple_assigns,2));

for t = 1:size(simple_assigns,2)
    
    for i = 1:length(waternum_R)
        
        water_com_assigns(watershed==waternum_R(i),t) = simple_assigns(length(waternum_L)+i,t);
    end
end

save(gifti(water_com_assigns),[waterdir '/avgcorrofcorr_allgrad_R_smooth255_edges_avg_uc_smooth_minima' num2str(minimadist) '_iter200_frach1_watershed_rawassn_minsize4.func.gii'])


%% Assign communities to minima rois, LR together
waterdir = '/data/cn4/laumannt/watershed_network/FCPROCESS_NEW/';

minima_rois = gifti([waterdir '/avgcorrofcorr_smooth255_allgrad_L_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center_rois.func.gii']);
minima_rois = minima_rois.cdata;
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
minima_rois(logical(mask),:) = 0;
roinum_L = size(minima_rois,2);

water_com_assigns = zeros(32492,size(simple_assigns,2));

for t = 1:size(simple_assigns,2)
    
    for i = 1:roinum_L
        
        water_com_assigns(logical(minima_rois(:,i)),t) = simple_assigns(i,t);
    end
end
save(gifti(water_com_assigns),'avgcorrofcorr_smooth255_allgrad_L_smooth255_edge_avg_uc_smooth_minima5_parcel_center_rois_rawassn_minsize4.func.gii')
                 

minima_rois = gifti([waterdir '/avgcorrofcorr_smooth255_allgrad_R_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center_rois.func.gii']);
minima_rois = minima_rois.cdata;
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
minima_rois(logical(mask)) = 0;
roinum_R = size(minima_rois,2);

water_com_assigns = zeros(32492,size(simple_assigns,2));

for t = 1:size(simple_assigns,2)
    
    for i = 1:roinum_R
        
        water_com_assigns(logical(minima_rois(:,i)),t) = simple_assigns(roinum_L+i,t);
    end
end

save(gifti(water_com_assigns),'avgcorrofcorr_smooth255_allgrad_R_smooth255_edge_avg_uc_smooth_minima5_parcel_center_rois_rawassn_minsize4.func.gii')


%% Assign communities to watersheds

water_com_assigns = zeros(32492,size(simple_assigns,2));

for t = 1:size(simple_assigns,2)
    
    for i = 1:length(waternum)
        
        water_com_assigns(watershed==waternum(i),t) = simple_assigns(i,t);
    end
end

%%

minima_com = zeros(32492,10);

for t = 1:size(minima_com,2)
    for n = 1:size(minima_ROI,2)
        
        minima_com(logical(minima_ROI(:,n)),t) = simple_assigns(n,t);
    end
end
    






