
waterdir = '/data/cn4/evan/RestingState/FC_Mapping_120/';
watershednames = {'120_L_wateredgethresh_watershedmerge_0.45.func.gii','120_R_wateredgethresh2_watershedmerge_0.45.func.gii'};%,'Poldrome_R_crossthresh_watershedmerge.func.gii'};
corrmatname = 'corrmat.mat';
roifilename = 'parcel_center.roi';
outputfolder = [waterdir '/120_LR_clxd_infomap/'];

kdenthresholds = [.005 .1];
kdeninterval = .0025;
distanceexclusion = 30;
networksizeminimum = 5;

outputstem = '120_';

HEMS = {'L','R'};
%hemname = {'LEFT';'RIGHT'};

calc_corrmat = 0;


medial_wallL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);   
    medial_wall{1} = ~(medial_wallL.cdata);  
    medial_wallR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);   
    medial_wall{2} = ~(medial_wallR.cdata);  
    
    %cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST_mod.txt';
     cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_DATALIST.txt';
 [subjects surfdatafile] = textread(cohortfile,'%s %s');
% %tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST_mod.txt';
 tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_TMASKLIST.txt';
 [subjects tmasks] = textread(tmaskfile,'%s %s');

 
 mkdir(outputfolder)
%timeseriesdataname = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_fieldmap_222/gradients_cifti_voxelclean_erode_concat_wateredge_correct/all34/allsubs_LR_timeseries.func.gii';
%allsubstmask = '/data/hcp-bluearc/home/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_fieldmap_222/gradients_cifti_voxelclean_erode_concat_wateredge_correct/all34/allsubs_total_tmask.txt';
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
         
         
        
    
    for hem = 1:length(HEMS)
        
        surf_timecourses = load(['/data/hcp-bluearc/home/laumannt/120_parcellation/surf_timecourses/' subjects{s} '_' HEMS{hem} '_time_dil10_32k_fs_LR_smooth2.55_noHEAD.func.gii']);

         surf_timecourses(:,1) = [];
% 
         surf_timecourses = surf_timecourses(:,logical(tmask));
         surf_timecourses(isnan(surf_timecourses)) = 0;
        
        
        
        thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']); 
        thismedialwall = ~thismedialwall.cdata;
        
        watershed = gifti([waterdir '/' watershednames{hem}]);
        watershed = watershed.cdata;
        
        waternum = unique(watershed);
        waternum(waternum==0) = [];
        
%         giftispace_surf_timecourses = zeros(32492,size(surf_timecourses,2));
%         
%         if strcmp(HEMS{hem},'L')
%             indicesincifti = 1:nnz(medial_wall{1});
%         else
%             indicesincifti = nnz(medial_wall{1})+1 : nnz(medial_wall{1})+nnz(medial_wall{2});
%         end
%         
%         giftispace_surf_timecourses(logical(thismedialwall),:) = surf_timecourses(indicesincifti,:);
            
giftispace_surf_timecourses = surf_timecourses; clear surf_timecourses

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

cd(waterdir)
save([outputfolder '/' corrmatname],'all_water_corrmat')

else
    load([outputfolder '/' corrmatname])
end


%% Find watershed parcel centers

center_coords_both = [];

xdistancemat = zeros(size(all_water_corrmat));

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
    
    
    
    %Save out midthickness coordinates of centers
    midthick = gifti([surfdir '/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii']);
    center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
    
    parcelcenterinds{hem} = indpos;
    center_coords_both = [center_coords_both ; center_coords];
    clear indpos
    
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' HEMS{hem} '.mat'])
    surf_distances{hem} = geo_distances;
    clear geo_distances
end

for hem = 1:length(HEMS)
    otherhem = abs(hem-3);
    
    heminds = (1:length(parcelcenterinds{hem})) + ((hem==2) * length(parcelcenterinds{1}));
    otherheminds = (1:length(parcelcenterinds{otherhem})) + ((hem==1) * length(parcelcenterinds{1}));
    
    
    hemxdistancemat = surf_distances{hem}(parcelcenterinds{hem},parcelcenterinds{hem});
    hemxdistancemat = hemxdistancemat > distanceexclusion;
    xdistancemat(heminds,heminds) = hemxdistancemat;
    
    
    hemxdistancemat = surf_distances{otherhem}(parcelcenterinds{hem},parcelcenterinds{otherhem});
    hemxdistancemat = hemxdistancemat > distanceexclusion;
    xdistancemat(heminds,otherheminds) = hemxdistancemat;
    
end
clear surf_distances
quickroifile(center_coords_both,[outputfolder '/' roifilename])


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
clear waternum


     ciftitemplatefile = ['/data/hcp-bluearc/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
    templatedata = gifti(ciftitemplatefile);
    ciftidata_towrite = zeros(length(templatedata.cdata),size(simple_assigns,2));


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
    
        ciftidata_towrite(((1:nnz(mask==0)) + (nnz(medial_wall{1})*(hem-1))) , t) = water_com_assigns((mask==0),t);
    
        
end
save(gifti(single(water_com_assigns)),[outputstem HEMS{hem} '_minsize' num2str(networksizeminimum) '.func.gii']);
end
      
 
     cifti_write_wHDR(ciftidata_towrite,ciftitemplatefile,[outputstem 'LR_minsize' num2str(networksizeminimum)])
 
    
    
 