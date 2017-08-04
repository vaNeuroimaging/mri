hem = 'L';
datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
[subjects subdata] = textread(datalist,'%s%s');

for s = 1:length(subjects)
    metricname = ['/data/cn4/evan/RestingState/Ind_variability/Subjects/' subjects{s} '/' subjects{s} '_L.func.gii'];

%metricnames = {'Temp.func.gii'};
% {'/data/cn4/evan/RestingState/FC_Mapping_120/120_L.func.gii',...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/C1/C1_L.func.gii',...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/C2/C2_L.func.gii',...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/C3/C3_L.func.gii',...
%     '/data/cn4/evan/RestingState/Ind_variability/Subjects/vc34120/vc34120_L.func.gii',...
%     '/data/cn4/evan/RestingState/Ind_variability/Subjects/vc33813/vc33813_L.func.gii',...
%     '/data/cn4/evan/RestingState/Ind_variability/Subjects/vc33769/vc33769_L.func.gii',...
%     '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_L.func.gii'};

midthicksurf = ['/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii'];
medialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.' hem '.32k_fs_LR.func.gii']); medialwall = ~medialwall.cdata;
%for metricnum = 1:length(metricnames)
    
evalc(['!wb_command -metric-gradient ' midthicksurf ' ' metricname ' Tempgrad.func.gii']);

metric = gifti(metricname); metric = metric.cdata(logical(medialwall));
gradmetric = gifti('Tempgrad.func.gii'); gradmetric = gradmetric.cdata(logical(medialwall));

smoothness = std(gradmetric) / std(metric);

disp([metricname ': ' num2str(smoothness)])
end