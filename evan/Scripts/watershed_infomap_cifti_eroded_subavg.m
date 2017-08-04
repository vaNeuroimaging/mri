
waterdir = '/data/cn4/evan/RestingState/FC_Mapping_120/';
watershednames = {'120_wateredge_L_crossthresh_watershedmerge.func.gii','120_wateredge_R_crossthresh_watershedmerge.func.gii'};
corrmatname = 'corrmat.mat';
roifilename = 'parcel_center.roi';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/120_wateredge_BOTH_crossthresh_infomap/';

kdenthresholds = [.005 .1];
kdeninterval = .005;
distanceexclusion = 20;
networksizeminimum = 4;

outputstem = '120_wateredge_BOTH_crossthresh';

HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};

medial_wallL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.L.32k_fs_LR.func.gii']);   
    medial_wall{1} = ~(medial_wallL.cdata);  
    medial_wallR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.R.32k_fs_LR.func.gii']);   
    medial_wall{2} = ~(medial_wallR.cdata);  
    
    %cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST_mod.txt';
    cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
[subjects surfdatafile] = textread(cohortfile,'%s %s');
%tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST_mod.txt';
tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

%% Watershed-based networks, Left and Right




for s = 1:length(subjects)
    clear watertime
    disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        watertime_both = [];
    
        timename = surfdatafile{s};%[subjects{s} '_cifti_timeseries_goodvoxels/vc32347_BOLD_LR_surf_subcort_smooth2.55_32k_fsLR'];
        evalc(['!wb_command -cifti-convert -to-gifti-ext ' timename ' ' waterdir '/Temp2.func.gii']);
        surf_timecourses = gifti([waterdir '/Temp2.func.gii']); surf_timecourses = surf_timecourses.cdata;
%         evalc(['!caret_command64 -file-convert -format-convert ASCII ' waterdir '/Temp2.func.gii']);
%         evalc(['!awk ''NF > 25'' ' waterdir '/Temp2.func.gii > ' waterdir '/Temp2_noHEAD.func.gii']);
        tmask = load(tmasks{s});
%         surf_timecourses = load([waterdir '/Temp2_noHEAD.func.gii']);
%         surf_timecourses(:,1) = [];

        surf_timecourses = surf_timecourses(:,logical(tmask));
        %nonnan_timepoints = ~logical(sum(isnan(surf_timecourses),1));
        %surf_timecourses = surf_timecourses(:,nonnan_timepoints);
        
    
    for hem = 1:length(watershednames)
        watershed = gifti([waterdir '/' watershednames{hem}]);
        watershed = watershed.cdata;
        
        waternum = unique(watershed);
        waternum(waternum==0) = [];
        
        giftispace_surf_timecourses = zeros(32492,size(surf_timecourses,2));
        
        if hem == 1
            indicesincifti = 1:nnz(medial_wall{1});
        else
            indicesincifti = nnz(medial_wall{1})+1 : nnz(medial_wall{1})+nnz(medial_wall{2});
        end
        
        giftispace_surf_timecourses(logical(medial_wall{hem}),:) = surf_timecourses(indicesincifti,:);
            
        for i = 1:length(waternum)
            
            waterind = find(watershed==waternum(i));
            watervert_timecourses = giftispace_surf_timecourses(waterind,:);
%            watervert_timecourses = watervert_timecourses(~logical(sum(isnan(watervert_timecourses),2)),:);
            watertime(:,i) = mean(watervert_timecourses,1)';
        end
        
        %delete([waterdir '/' timename '_noHEAD.func.gii'])
        %delete([waterdir '/' timename '.func.gii'])
        
        watertime_both = [watertime_both watertime];
        clear watertime
        
    end
    
    %watertime_both = [watertime{1} watertime{2}];
    %watertime_both = watertime{1};
    water_corrmat(:,:,s) = paircorr_mod(watertime_both);
    clear watertime_both

end
   


%all_water_corrmat = reshape(water_corrmat,[size(watertime_both,2) size(watertime_both,2) 120]);
all_water_corrmat = water_corrmat;

cd(waterdir)
save([waterdir '/' corrmatname],'all_water_corrmat')


%% Find watershed parcel centers

center_coords_both = [];

for hem = 1:length(watershednames)
    clear indpos
    surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    sphere = gifti([surfdir '/Conte69.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii']);
    
    [phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    
    
    water = gifti([waterdir '/' watershednames{hem}]);
    water = water.cdata;
    
    %Mask watersheds
%     maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
%     mask = gifti(maskname);
%     mask = mask.cdata;
    water(~logical(medial_wall{hem})) = 0;
    
    
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
%     maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
%     mask = gifti(maskname);
%     mask = mask.cdata;
    metric(~logical(medial_wall{hem})) = 0;   
    indpos = find(metric);
    
    %Save out midthickness coordinates of centers
    midthick = gifti([surfdir '/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    
    center_coords_both = [center_coords_both ; center_coords];
    %save(gifti(single(metric)),[dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.func.gii'])
    %dlmwrite([dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.txt'],indpos,'delimiter','\n')
    %save(['avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center.txt'],'indpos')
end

%center_coords_both = [center_coords{1};center_coords{2}];
%center_coords_both = center_coords{1};
quickroifile(center_coords_both,[waterdir '/' roifilename])

%% Make prmfile

prmfilename = [waterdir '/prmfile.txt'];
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,'watershed_bothhem','-append','delimiter','');
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
save(gifti(single(water_com_assigns)),[outputstem 'L_minsize' num2str(networksizeminimum) '.func.gii']);

      
if length(watershednames) > 1
    
    
    watershed = gifti([waterdir '/' watershednames{2}]);
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
            %water_com_assigns(watershed==waternum_R(i),t) = simple_assigns(i,t);
        end
    end
    
    save(gifti(single(water_com_assigns)),[outputstem 'R_minsize' num2str(networksizeminimum) '.func.gii']);
    
    gifti_to_cifti([outputstem 'L_minsize' num2str(networksizeminimum) '.func.gii'],[outputstem 'R_minsize' num2str(networksizeminimum) '.func.gii'],[outputstem 'BOTH_minsize' num2str(networksizeminimum)])
end