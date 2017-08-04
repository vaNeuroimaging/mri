
waterdir = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/';
watershednames = {'Poldrome_parcels_L.func.gii','Poldrome_parcels_R.func.gii'};
%{'120_L_unsmoothed_watershedmerge.func.gii','120_R_unsmoothed_watershedmerge.func.gii'};
corrmatname = 'corrmat.mat';
roifilename = 'parcel_center.roi';
outputfolder = [waterdir 'Poldrome_parcels_LR_infomap/'];

kdenthresholds = [.01 .05];
kdeninterval = .001;
distanceexclusion = 20;
networksizeminimum = 5;

outputstem = 'Poldrome_parcels_LR_assignments';


HEMS = {'L';'R'};
%hemname = {'LEFT';'RIGHT'};

calc_corrmat = 0;

medial_wallL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);   
    medial_wall{1} = ~(medial_wallL.cdata);  
    medial_wallR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);   
    medial_wall{2} = ~(medial_wallR.cdata);  
    
    %cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST_mod.txt';
%     cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
% [subjects surfdatafile] = textread(cohortfile,'%s %s');
% %tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST_mod.txt';
% tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
% [subjects tmasks] = textread(tmaskfile,'%s %s');

%timeseriesdataname = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_fieldmap_222/gradients_cifti_voxelclean_erode_concat_wateredge_correct/all34/allsubs_LR_timeseries.func.gii';
%allsubstmask = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_fieldmap_222/gradients_cifti_voxelclean_erode_concat_wateredge_correct/all34/allsubs_total_tmask.txt';
timeseriesdataname = '/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025_mirpad/cifti_correlation_normalwall/84sub_333_all_startpos50/poldrack_texas_startpos50.allsubs_LR_timeseries.dtseries.nii';
allsubstmask = '/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025_mirpad/cifti_correlation_normalwall/84sub_333_all_startpos50/allsubs_total_tmask.txt';
mkdir(outputfolder)
%% Watershed-based networks, Left and Right


if calc_corrmat

%for s = 1:length(subjects)
    %clear watertime
    %disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        watertime_both = [];
    
%         timename = surfdatafile{s};%[subjects{s} '_cifti_timeseries_goodvoxels/vc32347_BOLD_LR_surf_subcort_smooth2.55_32k_fsLR'];
%         evalc(['!wb_command -cifti-convert -to-gifti-ext ' timename ' ' waterdir '/Temp2.func.gii']);
%         surf_timecourses = gifti([waterdir '/Temp2.func.gii']); surf_timecourses = surf_timecourses.cdata;
% %         evalc(['!caret_command64 -file-convert -format-convert ASCII ' waterdir '/Temp2.func.gii']);
% %         evalc(['!awk ''NF > 25'' ' waterdir '/Temp2.func.gii > ' waterdir '/Temp2_noHEAD.func.gii']);
%         tmask = load(tmasks{s});
% %         surf_timecourses = load([waterdir '/Temp2_noHEAD.func.gii']);
% %         surf_timecourses(:,1) = [];
% 
%         surf_timecourses = surf_timecourses(:,logical(tmask));
        %nonnan_timepoints = ~logical(sum(isnan(surf_timecourses),1));
        %surf_timecourses = surf_timecourses(:,nonnan_timepoints);
        
        %surf_timecourses = gifti(timeseriesdataname); surf_timecourses = surf_timecourses.cdata;
        surf_timecourses = ft_read_cifti_mod(timeseriesdataname);
        surf_timecourses = surf_timecourses.data;
        tmask = load(allsubstmask);
        surf_timecourses = surf_timecourses(:,logical(tmask));
    
    for hem = 1:length(HEMS)
        
        thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']); 
        thismedialwall = ~thismedialwall.cdata;
        
        watershed = gifti([waterdir '/' watershednames{hem}]);
        watershed = watershed.cdata;
        
        waternum = unique(watershed);
        waternum(waternum==0) = [];
        
        giftispace_surf_timecourses = zeros(32492,size(surf_timecourses,2));
        
        if strcmp(HEMS{hem},'L')
            indicesincifti = 1:nnz(medial_wall{1});
        else
            indicesincifti = nnz(medial_wall{1})+1 : nnz(medial_wall{1})+nnz(medial_wall{2});
        end
        
        giftispace_surf_timecourses(logical(thismedialwall),:) = surf_timecourses(indicesincifti,:);
            
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
    water_corrmat = paircorr_mod(watertime_both);
    clear watertime_both

%end
   


%all_water_corrmat = reshape(water_corrmat,[size(watertime_both,2) size(watertime_both,2) 120]);
all_water_corrmat = water_corrmat;

cd(waterdir)
save([outputfolder '/' corrmatname],'all_water_corrmat')

else
   load([outputfolder '/' corrmatname]) 
end


%% Find watershed parcel centers

center_coords_both = [];

xdistancemat = ones(size(all_water_corrmat));

for hem = 1:length(HEMS)
    clear indpos
    surfdir = '/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/freesurfer_washu/FREESURFER_fs_LR/sub013/7112b_fs_LR/fsaverage_LR32k/';
    midthick = gifti([surfdir '/sub013.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    
    [phi theta r] = cart2sph(midthick.vertices(:,1), midthick.vertices(:,2),midthick.vertices(:,3));
    
    thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']); 
    thismedialwall = ~thismedialwall.cdata;
    
    water = gifti([waterdir '/' watershednames{hem}]);
    water = water.cdata;
    
    %Mask watersheds
%     maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
%     mask = gifti(maskname);
%     mask = mask.cdata;
    %water(~logical(thismedialwall)) = 0;
    
    
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
        
        
        meanX = mean(midthick.vertices(ind,1));
        meanY = mean(midthick.vertices(ind,2));
        meanZ = mean(midthick.vertices(ind,3));
        
        
        coord = [meanX meanY meanZ];
        midthick_coords = [midthick.vertices(ind,1) midthick.vertices(ind,2) midthick.vertices(ind,3)];
        
        rep_coord = repmat(coord, [size(midthick_coords,1) 1]);
        
        dist_coord = sum((midthick_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        indpos(w) = ind(indval);
        
    end
    
    metric = zeros(32492,1);    
    metric(indpos) = 1;
    
    %Mask metric
%     maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
%     mask = gifti(maskname);
%     mask = mask.cdata;
%     metric(~logical(thismedialwall)) = 0;   
%     indpos = find(metric);
    
    %Save out midthickness coordinates of centers
    midthick = gifti([surfdir '/sub013.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    
    load([surfdir 'geodesic_distance_' HEMS{hem} '.mat'])
    eval(['geo_distances = geodesic_distance_' HEMS{hem} '; clear geodesic_distance_' HEMS{hem}]);
    hemxdistancemat = geo_distances(indpos,indpos);
    clear geo_distances indpos
    hemxdistancemat = hemxdistancemat > distanceexclusion;
    
    if hem==1
        heminds = 1:size(hemxdistancemat,1); 
    else
        heminds = (size(xdistancemat,1) -size(hemxdistancemat,1) +1) : size(xdistancemat,1);
    end
    xdistancemat(heminds,heminds) = hemxdistancemat;
    
    center_coords_both = [center_coords_both ; center_coords];
    %save(gifti(single(metric)),[dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.func.gii'])
    %dlmwrite([dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.txt'],indpos,'delimiter','\n')
    %save(['avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center.txt'],'indpos')
end


%center_coords_both = [center_coords{1};center_coords{2}];
%center_coords_both = center_coords{1};
roifilename = [outputfolder '/' roifilename];
quickroifile(center_coords_both,roifilename)


%% Make prmfile

prmfilename = [outputfolder '/prmfile.txt'];
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,'watershed','-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
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

%% Modify color assignments with miminum network size criterion and regularize

filesinoutputfolder = dir(outputfolder);
for i=3:length(filesinoutputfolder);
    folderdatecreated(i-2) = filesinoutputfolder(i).datenum;
end

[maxval maxindex] = max(folderdatecreated);
trueoutputfolder = [outputfolder '/' filesinoutputfolder(maxindex+2).name];

cd(trueoutputfolder)

simple_assigns = modify_clrfile('simplify','rawassn.txt',networksizeminimum);
simple_assigns_reg = rawoutput2clr(simple_assigns);
simple_assigns_reg(simple_assigns_reg<=1) = 0;
simple_assigns_reg = simple_assigns_reg-1;
dlmwrite(['rawassn_minsize' num2str(networksizeminimum) '_regularized.txt'],simple_assigns_reg)

%% Assign communities to watersheds, LR together
clear waternum
simple_assigns_reg = load(['rawassn_minsize' num2str(networksizeminimum) '_regularized_consensus.txt']);

if length(HEMS) > 1
    ciftitemplatefile = ['/data/hcp-zfs/shared-nil/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
    ciftistruct = ft_read_cifti_mod(ciftitemplatefile); 
    %templatedata = gifti(ciftitemplatefile);
    ciftidata_towrite = zeros(size(ciftistruct.data),size(simple_assigns_reg,2));
    ciftistruct.data = [];
end

for hem = 1:length(HEMS)

watershed = gifti([waterdir '/' watershednames{hem}]);
watershed = watershed.cdata;
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
watershed(logical(mask)) = 0;
waternum{hem} = unique(watershed);
waternum{hem}(waternum{hem}==0) = [];

water_com_assigns = zeros(32492,size(simple_assigns_reg,2));

for t = 1:size(simple_assigns_reg,2)
    
    for i = 1:length(waternum{hem})
        
        water_com_assigns(watershed==waternum{hem}(i),t) = simple_assigns_reg(i + (length(waternum{1})*(hem==2)),t);
    end
    if length(HEMS) > 1
        ciftidata_towrite([1:nnz(mask==0)] + (nnz(medial_wall{1})*(hem-1)) , t) = water_com_assigns(mask==0,t);
    end
        
end
%save(gifti(single(water_com_assigns)),[outputstem HEMS{hem} '_minsize' num2str(networksizeminimum) '.func.gii']);
end
      
 if length(HEMS) > 1
     ciftistruct.data = ciftidata_towrite;
     ft_write_cifti_mod([outputstem '_minsize' num2str(networksizeminimum) '_regularized_consensus'],ciftistruct);
     %cifti_write_wHDR(ciftidata_towrite,ciftitemplatefile,[outputstem 'LR_minsize' num2str(networksizeminimum)])
 end
 
%  %% Tally-based consensus
%  tallymat_infomap(['rawassn_minsize' num2str(networksizeminimum) '.txt'],roifilename)
%  simple_assigns = load('tallymat_assn.txt');
%  
%  if length(HEMS) > 1
%     ciftitemplatefile = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
%     templatedata = gifti(ciftitemplatefile);
%     ciftidata_towrite = zeros(length(templatedata.cdata),size(simple_assigns,2));
% end
% 
% for hem = 1:length(HEMS)
% 
% watershed = gifti([waterdir '/' watershednames{hem}]);
% watershed = watershed.cdata;
% maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
% mask = gifti(maskname);
% mask = mask.cdata;
% watershed(logical(mask)) = 0;
% waternum{hem} = unique(watershed);
% waternum{hem}(waternum{hem}==0) = [];
% 
% water_com_assigns = zeros(32492,size(simple_assigns,2));
% 
% for t = 1:size(simple_assigns,2)
%     
%     for i = 1:length(waternum{hem})
%         
%         water_com_assigns(watershed==waternum{hem}(i),t) = simple_assigns(i + (length(waternum{1})*(hem==2)),t);
%     end
%     if length(HEMS) > 1
%         ciftidata_towrite([1:nnz(mask==0)] + (nnz(medial_wall{1})*(hem-1)) , t) = water_com_assigns(mask==0,t);
%     end
%         
% end
% %save(gifti(single(water_com_assigns)),[outputstem HEMS{hem} '_minsize' num2str(networksizeminimum) '_tallymatconsensus.func.gii']);
% end
%       
%  if length(HEMS) > 1
%      cifti_write_wHDR(ciftidata_towrite,ciftitemplatefile,[outputstem 'LR_minsize' num2str(networksizeminimum) '_tallymatconsensus'])
%  end
%  