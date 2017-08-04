xdistance = 30;

kden_int = .002;
kdens = [.006 : kden_int : .02];

corrmatname = '/data/hcp-bluearc/home/laumannt/120_parcellation/modified_cifti_network/cifti_normalwall/BOTH/cifti_avgcrosscorr.mat';
roifilename = '/data/hcp-bluearc/home/laumannt/120_parcellation/modified_cifti_network/cifti_normalwall/BOTH/cifti_coords_LR.roi';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/120_vertexwise_homotopicthreshx/';

hems = {'L','R'};

%%
disp('Loading distance files')
for hemnum = 1:length(hems)
    hem = hems{hemnum};

thismask = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']); mask{hemnum} = ~thismask.cdata;
load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])
surf_distances{hemnum} = geo_distances;
clear geo_distances
end

subcort_verts = textread('/data/cn4/evan/ROIs/Subcort.txt');

allverts = load('/data/hcp-bluearc/home/laumannt/120_parcellation/modified_cifti_network/cifti_normalwall/BOTH/cifti_coords_LR.txt');
%%
disp('Calculating distance exclusion matrix')

xdistancemat = logical(ones(size(allverts,1)));

for hemnum = 1:length(hems)
    otherhemnum = abs(hemnum-3);
    
    heminds = (1:nnz(mask{hemnum})) + ((hemnum==2) * nnz(mask{1}));
    otherheminds = (1:nnz(mask{otherhemnum})) + ((hemnum==1) * nnz(mask{1}));
     
    
    hemxdistancemat = surf_distances{hemnum}(logical(mask{hemnum}),logical(mask{hemnum}));
    hemxdistancemat = hemxdistancemat > xdistance;
    xdistancemat(heminds,heminds) = logical(hemxdistancemat);
    
    
%      hemxdistancemat = surf_distances{otherhemnum}(logical(mask{hemnum}),logical(mask{otherhemnum}));
%      hemxdistancemat = hemxdistancemat > xdistance;
%      xdistancemat(heminds,otherheminds) = logical(hemxdistancemat);
    
     clear hemxdistancemat
end

for subcort_index = (size(allverts,1) - size(subcort_verts,1) + 1) : size(allverts,1)
    dist_fromvert = distance(allverts(subcort_index,:)',allverts');
    xdistancemat(subcort_index,:) = logical(dist_fromvert > xdistance);
    xdistancemat(:,subcort_index) = logical(dist_fromvert > xdistance);
end

%%
%disp('Loading corrmat file')
%load(corrmatname);
% avgcrosscorr_LR = avgcrosscorr_LR .* xdistancemat;
% clear xdistancemat surf_distances
% for kdennum = 1:length(kdens);
%     thiskden = kdens(kdennum);
%     [rmat rs(kdennum) kden] = matrix_thresholder_onlyrs(avgcrosscorr_LR,thiskden,'kden');
%     
% end
% 
% rint_approx = mean(diff(rs));

% %%
% disp('Reloading distance files and recalculating distance exclusion matrix')
% for hemnum = 1:length(hems)
%     hem = hems{hemnum};
% 
% load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])
% surf_distances{hemnum} = geo_distances;
% clear geo_distances
% end
% 
% xdistancemat = logical(ones(size(allverts,1)));
% 
% for hemnum = 1:length(hems)
%     otherhemnum = abs(hemnum-3);
%     
%     heminds = (1:nnz(mask{hemnum})) + ((hemnum==2) * nnz(mask{1}));
%     otherheminds = (1:nnz(mask{otherhemnum})) + ((hemnum==1) * nnz(mask{1}));
%      
%     
%     hemxdistancemat = surf_distances{hemnum}(logical(mask{hemnum}),logical(mask{hemnum}));
%     hemxdistancemat = hemxdistancemat > xdistance;
%     xdistancemat(heminds,heminds) = logical(hemxdistancemat);
%     
%      clear hemxdistancemat
% end
% for subcort_index = (size(allverts,1) - size(subcort_verts,1) + 1) : size(allverts,1)
%     dist_fromvert = distance(allverts(subcort_index,:)',allverts');
%     xdistancemat(subcort_index,:) = logical(dist_fromvert > xdistance);
%     xdistancemat(:,subcort_index) = logical(dist_fromvert > xdistance);
% end
% clear surf_distances



%% Make prmfile
disp('Running infomap')
mkdir(outputfolder)
prmfilename = [outputfolder '/prmfile.txt'];
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,'120_LR','-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdens(1)),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kden_int),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdens(end)),'-append','delimiter','');
dlmwrite(prmfilename,outputfolder,'-append','delimiter','');
dlmwrite(prmfilename,num2str(.1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(12),'-append','delimiter','');
dlmwrite(prmfilename,num2str(5),'-append','delimiter','');
dlmwrite(prmfilename,'kden','-append','delimiter','');

%% Run infomap
cd(outputfolder)
graphcluster_Evan_surfacexd([outputfolder '/prmfile.txt'],'thr',xdistance,xdistancemat,1,0,'infomap');


    