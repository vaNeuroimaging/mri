thisdir = pwd;
subjectnums = [];

kdenthresholds = [.005 .1];
kdeninterval = .005;
distanceexclusion = 20;
networksizeminimum = 3;

outputstem = 'Parcels_';

HEMS = {'L'};%;'R'};
%hemname = {'LEFT';'RIGHT'};

medial_wallL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.L.32k_fs_LR.func.gii']);   
    medial_wall{1} = ~(medial_wallL.cdata);  
    medial_wallR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.R.32k_fs_LR.func.gii']);   
    medial_wall{2} = ~(medial_wallR.cdata);  
    
 
datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';
    
    [subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

if ~isempty(subjectnums)
    subjects = subjects(subjectnums); subdata = subdata(subjectnums); tmasks = tmasks(subjectnums);
end

corrmatname = 'corrmat.mat';
roifilename = 'parcel_center.roi';

for s = 1:length(subjects)
    
    

%waterdir = ['/data/cn4/evan/RestingState/FC_Mapping_120/'];
%watershednames = {['120_L_target240_watershedmerge_noratio.func.gii']};%,'Poldrome_R_crossthresh_watershedmerge.func.gii'};
waterdir = ['/data/cn4/evan/RestingState/Ind_variability/Subjects/' subjects{s} '/'];
watershednames = {[subjects{s} '_L_target240_watershedmerge_noratio.func.gii']};%,'Poldrome_R_crossthresh_watershedmerge.func.gii'};


%outputfolder = ['/data/cn4/evan/RestingState/Ind_variability/Subjects/' subjects{s} '/Parcels_group_target240_infomap/'];
outputfolder = [waterdir '/Parcels_target240_noratio_infomap/'];
evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{s} ' ' thisdir '/Temp1.func.gii']);
timeseriesdataname = [thisdir '/Temp1.func.gii'];
allsubstmask = tmasks{s};

%% Watershed-based networks, Left and Right




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
        
        surf_timecourses = gifti(timeseriesdataname); surf_timecourses = surf_timecourses.cdata;
        tmask = load(allsubstmask);
        surf_timecourses = surf_timecourses(:,logical(tmask));
    
    for hem = 1:length(HEMS)
        
        signalstrength = gifti(['/data/cn4/evan/ROIs/all_meanimage_' HEMS{hem} '.func.gii']);
        signalstrength = signalstrength.cdata;
        
        thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.' HEMS{hem} '.32k_fs_LR.func.gii']); 
        thismedialwall = ~thismedialwall.cdata;
        
        watershed = gifti([waterdir '/' watershednames{hem}]);
        watershed = watershed.cdata;
        watershed(signalstrength<750) = 0;
        
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
mkdir(outputfolder)
cd(outputfolder)
save([outputfolder '/' corrmatname],'all_water_corrmat')


%% Find watershed parcel centers

center_coords_both = [];

for hem = 1:length(HEMS)
    clear indpos
    surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    sphere = gifti([surfdir '/Conte69.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii']);
    signalstrength = gifti(['/data/cn4/evan/ROIs/all_meanimage_' HEMS{hem} '.func.gii']);
        signalstrength = signalstrength.cdata;
    
    [phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));
    
    thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.' HEMS{hem} '.32k_fs_LR.func.gii']); 
    thismedialwall = ~thismedialwall.cdata;
    
    water = gifti([waterdir '/' watershednames{hem}]);
    water = water.cdata;
    water(signalstrength<750) = 0;
    
    %Mask watersheds
%     maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
%     mask = gifti(maskname);
%     mask = mask.cdata;
    water(~logical(thismedialwall)) = 0;
    
    
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
%     metric(~logical(thismedialwall)) = 0;   
%     indpos = find(metric);
    
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
quickroifile(center_coords_both,[outputfolder '/' roifilename])

%% Make prmfile

prmfilename = [outputfolder '/prmfile.txt'];
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,'watershed_bothhem','-append','delimiter','');
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

graphcluster_Evan([outputfolder '/prmfile.txt'],'thr',distanceexclusion  ,1,0,'infomap');

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

for hem = 1:length(HEMS)
    signalstrength = gifti(['/data/cn4/evan/ROIs/all_meanimage_' HEMS{hem} '.func.gii']);
        signalstrength = signalstrength.cdata;

watershed = gifti([waterdir '/' watershednames{hem}]);
watershed = watershed.cdata;
watershed(signalstrength<750) = 0;
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
watershed(logical(mask)) = 0;
parcelIDs{hem} = unique(watershed);
parcelIDs{hem}(parcelIDs{hem}==0) = [];

water_com_assigns = zeros(32492,size(simple_assigns,2));

for t = 1:size(simple_assigns,2)
    
    for i = 1:length(parcelIDs{hem})
        
        water_com_assigns(watershed==parcelIDs{hem}(i),t) = simple_assigns(i + (length(parcelIDs{1})*(hem==2)),t);
    end
end
save(gifti(single(water_com_assigns)),[outputstem HEMS{hem} '_minsize' num2str(networksizeminimum) '.func.gii']);
end
      
% if length(HEMS) > 1
%     gifti_to_cifti([outputstem 'L_minsize' num2str(networksizeminimum) '.func.gii'],[outputstem 'R_minsize' num2str(networksizeminimum) '.func.gii'],[outputstem 'BOTH_minsize' num2str(networksizeminimum)])
% end

end