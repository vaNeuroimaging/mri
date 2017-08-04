%%
waterdir = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/';
watershednames = {'Poldrome_subsurf_edge_L_watershedmerge_0.45.func.gii'};
%{'120_L_unsmoothed_watershedmerge.func.gii','120_R_unsmoothed_watershedmerge.func.gii'};
corrmatname = 'corrmat.mat';
roifilename = 'parcel_center.roi';
outputfolder = [waterdir 'Poldrome_subsurf_edge_L_watershedmerge_045_infomap/'];

kdenthresholds = [.01 .06];
kdeninterval = .001;
distanceexclusion = 30;
networksizeminimum = 2;

outputstem = 'Poldrome_subsurf_edge_L_watershedmerge_045_watershedassignments';





HEMS = {'L'};%,'R'};
%hemname = {'LEFT';'RIGHT'};

%alpha = .05;

calc_corrmat = 0;


medial_wallL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);   
    medial_wall{1} = ~(medial_wallL.cdata);  
    medial_wallR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);   
    medial_wall{2} = ~(medial_wallR.cdata);  
    

 tmaskfile = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/poldrome_allses_selected_final_TMASKLIST.txt';
 [subjects tmasks] = textread(tmaskfile,'%s %s');

 
 mkdir(outputfolder)
%timeseriesdataname = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_fieldmap_222/gradients_cifti_voxelclean_erode_concat_wateredge_correct/all34/allsubs_LR_timeseries.func.gii';
%allsubstmask = '/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_fieldmap_222/gradients_cifti_voxelclean_erode_concat_wateredge_correct/all34/allsubs_total_tmask.txt';
%timeseriesdataname = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_LR_timeseries.func.gii';
%allsubstmask = '/data/cn4/evan/RestingState/FC_Mapping_120/allsubs_total_tmask.txt';

%% Watershed-based networks, Left and Right

if calc_corrmat


for s = 1:length(subjects)
    %clear watertime
    disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        watertime_both = [];
    
%         timename = surfdatafile{s};%[subjects{s} '_cifti_timeseries_goodvoxels/vc32347_BOLD_LR_surf_subcort_smooth2.55_32k_fsLR'];
%         evalc(['!wb_command -cifti-convert -to-gifti-ext ' timename ' ' waterdir '/Temp2.func.gii']);
%         surf_timecourses = gifti([waterdir '/Temp2.func.gii']); surf_timecourses = surf_timecourses.cdata;
% %         evalc(['!caret_command64 -file-convert -format-convert ASCII ' waterdir '/Temp2.func.gii']);
% %         evalc(['!awk ''NF > 25'' ' waterdir '/Temp2.func.gii > ' waterdir '/Temp2_noHEAD.func.gii']);
         tmask = load(tmasks{s});
         
        cifti_timeseries = cifti_read(['/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025/cifti_timeseries_normalwall/' subjects{s} '_BOLD_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
        cifti_timeseries = cifti_timeseries(:,logical(tmask));
        cifti_timeseries(isnan(cifti_timeseries)) = 0;
    
    for hem = 1:length(HEMS)
        
        %surf_timecourses = load(['/data/hcp-zfs/home/laumannt/120_parcellation/surf_timecourses/' subjects{s} '_' HEMS{hem} '_time_dil10_32k_fs_LR_smooth2.55_noHEAD.func.gii']);
        %/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_lininterp_222_sub018_reg_FD025/surf_timecourses/sub104_L_time_333_dil10_32k_fs_LR_smooth2.55.func.gii
        
        thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']); 
        thismedialwall = ~thismedialwall.cdata;
        
        watershed = gifti([waterdir '/' watershednames{hem}]);
        watershed = watershed.cdata;
        
        waternum = unique(watershed);
        waternum(waternum==0) = [];
        
        giftispace_surf_timecourses = zeros(32492,size(cifti_timeseries,2));
        giftispace_surf_timecourses(logical(thismedialwall),:) = cifti_timeseries([1:nnz(thismedialwall)] + (nnz(medial_wall{1}) * strcmp(HEMS{hem},'R')),:);


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
all_water_corrmat = mean(water_corrmat,3);
std_water_corrmat = std(water_corrmat,[],3);

cd(waterdir)
save([outputfolder '/' corrmatname],'all_water_corrmat')

cd(waterdir)
save([outputfolder '/std_' corrmatname],'std_water_corrmat')

else
    load([outputfolder '/' corrmatname])
    load([outputfolder '/std_' corrmatname])
end


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

dlmwrite(prmfilename,[outputfolder '/' corrmatname],'-append','delimiter','');
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


cd(outputfolder)
graphcluster_Evan_surfacexd([outputfolder '/prmfile.txt'],'thr',distanceexclusion,xdistancemat,1,0,'infomap');


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