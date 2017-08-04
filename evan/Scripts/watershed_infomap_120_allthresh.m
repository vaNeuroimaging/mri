
waterdir = '/data/cn5/selfRegulation/V4Process_nosmooth/gradients_120_108_combined_subsurf_nosmooth/';
%'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/';
%
watershednames = {'120_108_combined_L_watershedmerge_0.35_tweaked.func.gii','120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'};%,'Poldrome_R_crossthresh_watershedmerge.func.gii'};
%{'120_subsurf_L_nosmooth_watershedmerge_0.4_tweaked.func.gii','120_subsurf_R_nosmooth_watershedmerge_0.4_tweaked.func.gii'};
%
corrmatname = 'corrmat.mat';
roifilename = 'parcel_center.roi';
outputfolder = [waterdir '/120_108_combined_LR_infomap_in120/'];

kdenthresholds = [.0001 .05];

distanceexclusion = 30;
networksizeminimum = 4;

outputstem = '120_108_';

HEMS = {'L','R'};
%hemname = {'LEFT';'RIGHT'};

calc_corrmat = 0;
calc_roifile = 0;


medial_wallL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);   
    medial_wall{1} = ~(medial_wallL.cdata);  
    medial_wallR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);   
    medial_wall{2} = ~(medial_wallR.cdata);  
    
    %cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Surfdatalist_120.txt';
    cohortfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
 [subjects surfdatafile] = textread(cohortfile,'%s %s');
% %tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST_mod.txt';
 tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
 [subjects tmasks] = textread(tmaskfile,'%s %s');

 kdeninterval = 1; %meaningless
 
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
         
         
        
    
    
        
        %surf_timecourses = load(['/data/hcp-zfs/home/laumannt/120_parcellation/surf_timecourses/' subjects{s} '_' HEMS{hem} '_time_dil10_32k_fs_LR_smooth2.55_noHEAD.func.gii']);

         %surf_timecourses(:,1) = [];
% 
        surf_timecourses = cifti_read(surfdatafile{s});


         surf_timecourses = surf_timecourses(:,logical(tmask));
         surf_timecourses(isnan(surf_timecourses)) = 0;
        
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
            
%giftispace_surf_timecourses = surf_timecourses; clear surf_timecourses

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
    water_corrmat(:,:,s) = FisherTransform(paircorr_mod(watertime_both));
    clear watertime_both

end
   


%all_water_corrmat = reshape(water_corrmat,[size(watertime_both,2) size(watertime_both,2) 120]);
all_water_corrmat = mean(water_corrmat,3);

cd(waterdir)
save([outputfolder '/' corrmatname],'all_water_corrmat')

else
    load([outputfolder '/' corrmatname])
end


%% Find watershed parcel centers

if calc_roifile

center_coords_both = [];

xdistancemat = ones(size(all_water_corrmat));

for hem = 1:length(HEMS)
    clear indpos
    surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    midthick = gifti([surfdir '/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    
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
    midthick = gifti([surfdir '/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    
    if distanceexclusion > 0
        load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' HEMS{hem} '.mat'])
        hemxdistancemat = geo_distances(indpos,indpos);
        clear geo_distances indpos
        hemxdistancemat = hemxdistancemat > distanceexclusion;
        
        if hem==1
            heminds = 1:size(hemxdistancemat,1);
        else
            heminds = (size(xdistancemat,1) -size(hemxdistancemat,1) +1) : size(xdistancemat,1);
        end
        xdistancemat(heminds,heminds) = hemxdistancemat;
    end
    
    center_coords_both = [center_coords_both ; center_coords];
    %save(gifti(single(metric)),[dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.func.gii'])
    %dlmwrite([dir '/Edges_' HEMS{hem} '_minima3watermerged_parcel_center.txt'],indpos,'delimiter','\n')
    %save(['avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed_parcel_center.txt'],'indpos')
end


%center_coords_both = [center_coords{1};center_coords{2}];
%center_coords_both = center_coords{1};
quickroifile(center_coords_both,[outputfolder '/' roifilename])
save('xdistancemat.mat','xdistancemat')

else
    load([outputfolder '/xdistancemat.mat'])
end

%%
connectionvals = all_water_corrmat .* xdistancemat .* triu(ones(size(all_water_corrmat)),1);
connectionvals = connectionvals(:); connectionvals(connectionvals==0) = [];

sortedvals = sort(connectionvals,'descend');



allthresholds = unique(sortedvals(round(length(sortedvals).*kdenthresholds(1)) : round(length(sortedvals).*kdenthresholds(2)))) + .0000000001;
allthresholds = allthresholds([length(allthresholds) : -1 : 1]);

%% Make prmfile

prmfilename = [outputfolder '/prmfile.txt'];
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,'parcel_allthresh','-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(1)),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdeninterval),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(2)),'-append','delimiter','');
dlmwrite(prmfilename,outputfolder,'-append','delimiter','');
dlmwrite(prmfilename,num2str(.1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(12),'-append','delimiter','');
dlmwrite(prmfilename,num2str(5),'-append','delimiter','');
dlmwrite(prmfilename,'r','-append','delimiter','');

%% Run infomap
cd(outputfolder)
graphcluster_Evan_surfacexd_customthresholds([outputfolder '/prmfile.txt'],'thr',distanceexclusion,xdistancemat,0,0,'infomap',allthresholds);

%% Modify color assignments with miminum network size criterion

filesinoutputfolder = dir(outputfolder);
for i=3:length(filesinoutputfolder);
    folderdatecreated(i-2) = filesinoutputfolder(i).datenum;
end

[maxval maxindex] = max(folderdatecreated);
trueoutputfolder = [outputfolder '/' filesinoutputfolder(maxindex+2).name];

cd(trueoutputfolder)

simple_assigns = modify_clrfile('simplify','rawassn.txt',networksizeminimum);
%%
% %% Assign communities to watersheds, LR together
% clear waternum
% 
% simple_assigns = load(['rawassn_minsize' num2str(networksizeminimum) '.txt']);
% 
% if length(HEMS) > 1
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
% save(gifti(single(water_com_assigns)),[outputstem HEMS{hem} '_minsize' num2str(networksizeminimum) '.func.gii']);
% end
%       
%  if length(HEMS) > 1
%      cifti_write_wHDR(ciftidata_towrite,ciftitemplatefile,[outputstem 'LR_minsize' num2str(networksizeminimum)])
%  end
 
 %%
 if length(HEMS) > 1
    ciftitemplatefile = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
    templatedata = gifti(ciftitemplatefile);
    ciftidata_towrite = zeros(length(templatedata.cdata),1);
 end
 assns = load(['rawassn_minsize' num2str(networksizeminimum) '.txt']);
 tallymat_infomap(assns,[outputfolder '/' roifilename]);
simple_assigns = modify_clrfile_nopic('simplify','tallymat_assn.txt',2);
 

 for hem = 1:length(HEMS)

watershed = gifti([waterdir '/' watershednames{hem}]);
watershed = watershed.cdata;
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
watershed(logical(mask)) = 0;
waternum{hem} = unique(watershed);
waternum{hem}(waternum{hem}==0) = [];

water_com_assigns = zeros(32492,size(simple_assigns,2));

for t = 1:size(simple_assigns,2)
    
    for i = 1:length(waternum{hem})
        
        water_com_assigns(watershed==waternum{hem}(i),t) = simple_assigns(i + (length(waternum{1})*(hem==2)),t);
    end
    if length(HEMS) > 1
        ciftidata_towrite([1:nnz(mask==0)] + (nnz(medial_wall{1})*(hem-1)) , t) = water_com_assigns(mask==0,t);
    end
        
end
%save(gifti(single(water_com_assigns)),[outputstem HEMS{hem} '_minsize' num2str(networksizeminimum) '.func.gii']);
end
      
 if length(HEMS) > 1
     cifti_write_wHDR(ciftidata_towrite,ciftitemplatefile,['tallyconsensus_' outputstem 'LR_minsize2'])
 end
 
    
    
 