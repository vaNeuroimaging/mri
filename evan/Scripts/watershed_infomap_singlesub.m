%%
waterdir = '/data/cn4/evan/Temp/';
watershednames = {'AdjustedWatershedEdges_1.2_watershedmerge.func.gii'};
%{'120_L_unsmoothed_watershedmerge.func.gii','120_R_unsmoothed_watershedmerge.func.gii'};
corrmatname = 'corrmat.mat';
roifilename = 'parcel_center.roi';
outputfolder = [waterdir 'vc33416_L_watershedinfomap_thresh1p2/'];

kdenthresholds = [.015 .1];
kdeninterval = .005;
distanceexclusion = 20;
networksizeminimum = 2;

outputstem = 'vc33416_L_watershedinfomap_thresh1p2_assignments';

%% Watershed-based networks, Left and Right
cohs = {'C1';'C2';'C3'};
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};



surftimedir = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort/';

%cd(surftimedir)
    
% cohortdir = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients/120cohort'];
% [subjects tmasks] = textread([cohortdir '/NEW_nokids_TMASKLIST.txt'],'%q %q');

%[subjects tmasks] = textread('/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS/COHORTSELECT/NEW_TMASKLIST.txt','%s%s');
subjects{1} = 'vc33416';
tmasks{1} = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/FCPROCESS_bandpass_interp_nosmooth/vc33416/total_tmask.txt';

for s = 1:length(subjects)
    clear watertime
    disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        
    for hem = 1%:2
        watershed = gifti([waterdir '/' watershednames{hem}]);
        watershed = watershed.cdata;
        maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
        mask = gifti(maskname);
        mask = mask.cdata;
        watershed(logical(mask)) = 0;
        waternum = unique(watershed);
        waternum(waternum==0) = [];
        
        tic
        %timename = [surftimedir subjects{s} '_BOLD_' HEMS{hem} '_smooth2.55_32k_fsLR'];
        timename = [subjects{s} '_BOLD_' HEMS{hem} '_smooth2.55_32k_fsLR'];
        copyfile([surftimedir '/' timename '.func.gii'],[waterdir '/' timename '.func.gii']);
        evalc(['!caret_command64 -file-convert -format-convert ASCII ' waterdir '/' timename '.func.gii']);
        evalc(['!awk ''NF > 25'' ' waterdir '/' timename '.func.gii > ' waterdir '/' timename '_noHEAD.func.gii']);
        tmask = load(tmasks{s});
        surf_timecourses = load([waterdir '/' timename '_noHEAD.func.gii']);
        surf_timecourses(:,1) = [];
        surf_timecourses = surf_timecourses(:,logical(tmask));
        
        for i = 1:length(waternum)
            
            waterind = find(watershed==waternum(i));
            watertime{hem}(:,i) = mean(surf_timecourses(waterind,:))';
        end
        
        delete([waterdir '/' timename '_noHEAD.func.gii'])
        delete([waterdir '/' timename '.func.gii'])
        
    end
    
    watertime_both = [];
    for i = 1:length(watertime)
        watertime_both = [watertime_both watertime{i}];
    end
    water_corrmat(:,:,s) = paircorr_mod(watertime_both);

end
   


%all_water_corrmat = reshape(water_corrmat,[size(watertime_both,2) size(watertime_both,2) 120]);
all_water_corrmat = water_corrmat;

cd(waterdir)
save([waterdir '/' corrmatname],'all_water_corrmat')


%% Find watershed parcel centers

for hem = 1%:2
    clear indpos
    surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    sphere = gifti([surfdir '/Conte69.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii']);
    
    [phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    
    
    water = gifti([waterdir '/' watershednames{hem}]);
    water = water.cdata;
    
    %Mask watersheds
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = mask.cdata;
    water(logical(mask)) = 0;
    
    
    waternum = unique(water);
    waternum(waternum==0) = [];
    
    for w = 1:length(waternum)
        
        ind = find(water==waternum(w));
        
        %     meanphi = mean(phi(ind));
        %     meantheta = mean(theta(ind));
        %
        %     coord = [meanphi meantheta];
        %     sphere_coords = [phi(ind) theta(ind)];
        %
        %     rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        %
        %     dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        
        
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        
        
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
        
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        indpos(w) = ind(indval);
        
    end
    
    metric = zeros(32492,1);    
    metric(indpos) = 1;
    
    %Mask metric
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = mask.cdata;
    metric(logical(mask)) = 0;   
    indpos = find(metric);
    
    %Save out midthickness coordinates of centers
    midthick = gifti([surfdir '/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    center_coords{hem} = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    %save(gifti(single(metric)),[dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.func.gii'])
    %dlmwrite([dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.txt'],indpos,'delimiter','\n')
    %save(['avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center.txt'],'indpos')
end

center_coords_both = [];
for i = 1:length(center_coords)
    center_coords_both = [center_coords_both; center_coords{i}];
end
quickroifile(center_coords_both,[waterdir '/' roifilename])

%% Make prmfile

prmfilename = [waterdir '/prmfile.txt'];
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,'watershed','-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(length(subjects)),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(1)),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdeninterval),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(2)),'-append','delimiter','');
dlmwrite(prmfilename,outputfolder,'-append','delimiter','');
dlmwrite(prmfilename,num2str(.1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(12),'-append','delimiter','');
dlmwrite(prmfilename,num2str(5),'-append','delimiter','');
dlmwrite(prmfilename,'kden','-append','delimiter','');

%% Run infomap

graphcluster_Evan([waterdir '/prmfile.txt'],'thr',distanceexclusion  ,1,0,'infomap');

%% Modify color assignments with miminum network size criterion

filesinoutputfolder = dir(outputfolder);
for i=3:length(filesinoutputfolder);
    folderdatecreated(i-2) = filesinoutputfolder(i).datenum;
end

[maxval maxindex] = max(folderdatecreated);
trueoutputfolder = [outputfolder '/' filesinoutputfolder(maxindex+2).name];

cd(trueoutputfolder)

simple_assigns = modify_clrfile('simplify','rawassn.txt',networksizeminimum);

%% Assign communities to watersheds, LR together

watershed = gifti([waterdir '/' watershednames{1}]);
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
save(gifti(water_com_assigns),[outputstem 'L_minsize' num2str(networksizeminimum) '.func.gii']);

                 

% watershed = gifti([waterdir '/' watershednames{2}]);
% watershed = watershed.cdata;
% maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii'];
% mask = gifti(maskname);
% mask = mask.cdata;
% watershed(logical(mask)) = 0;
% waternum_R = unique(watershed);
% waternum_R(waternum_R==0) = [];
% 
% water_com_assigns = zeros(32492,size(simple_assigns,2));
% 
% for t = 1:size(simple_assigns,2)
%     
%     for i = 1:length(waternum_R)
%         
%         water_com_assigns(watershed==waternum_R(i),t) = simple_assigns(length(waternum_L)+i,t);
%         %water_com_assigns(watershed==waternum_R(i),t) = simple_assigns(i,t);
%     end
% end
% 
% save(gifti(water_com_assigns),[outputstem 'R_minsize' num2str(networksizeminimum) '.func.gii']);
% 
%  gifti_to_cifti([outputstem 'L_minsize' num2str(networksizeminimum) '.func.gii'],[outputstem 'R_minsize' num2str(networksizeminimum) '.func.gii'],[outputstem 'BOTH_minsize' num2str(networksizeminimum) '.func.gii'])