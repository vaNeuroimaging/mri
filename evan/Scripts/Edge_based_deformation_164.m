function Edge_based_deformation_164(sourceedgemap,targetedgemap,parcels,hem)

system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample ' sourceedgemap ' /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii ADAP_BARY_AREA ' sourceedgemap(1:end-9) '_164.func.gii'])
system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample ' targetedgemap ' /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii ADAP_BARY_AREA ' targetedgemap(1:end-9) '_164.func.gii'])
system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample ' parcels ' /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii ADAP_BARY_AREA ' parcels(1:end-9) '_164.func.gii'])

sourceedgemap = [sourceedgemap(1:end-9) '_164.func.gii'];
targetedgemap = [targetedgemap(1:end-9) '_164.func.gii'];
parcels = [parcels(1:end-9) '_164.func.gii'];

thisfolder = pwd;

slashloc = strfind(sourceedgemap,'/');
if isempty(slashloc)
    outfolder = pwd;
else
    outfolder = sourceedgemap(1:slashloc(end));
end;

cd(outfolder)

%hem = 'L';

if strcmp(hem,'L')
hemlow = 'l';
fullhem = 'left';
elseif strcmp(hem,'R')
hemlow = 'r';
fullhem = 'right';
end

% sourceedgemap = '';
% targetedgemap = '';

mkdir([outfolder '/source/'])
mkdir([outfolder '/target/'])

coordfile = ['/data/cn4/gagan/PARCELLATION_2012/TOOLS_N_STUFF/FSAVG2FSLR_SCRIPTS/global/templates/standard_mesh_atlases/fsaverage.' hem '_LR.spherical_std.164k_fs_LR.coord'];
topofile = ['/data/cn4/gagan/PARCELLATION_2012/TOOLS_N_STUFF/FSAVG2FSLR_SCRIPTS/global/templates/standard_mesh_atlases/fsaverage.' hem '.closed.164k_fs_LR.topo'];
sulcfile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.sulc.32k_fs_LR_' hem '.func.gii'];
whitefile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.white.32k_fs_LR.surf.gii'];

system(['caret_command64 -file-convert -sc -is CARET ' coordfile ' ' topofile ' -os FSS source/' hemlow 'h.sphere.asc'])
system(['mris_convert source/' hemlow 'h.sphere.asc source/' hemlow 'h.sphere'])
copyfile(['source/' hemlow 'h.sphere'],['target/' hemlow 'h.sphere'])

% system(['caret_command64 -file-convert -cs2fsc ' sulcfile ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.coord.gii source/' hemlow 'h.sulc.asc']);
% system(['mris_convert -c source/' hemlow 'h.sulc.asc source/' hemlow 'h.sphere.asc source/' hemlow 'h.sulc'])
% copyfile(['source/' hemlow 'h.sulc'],['target/' hemlow 'h.sulc'])


% system(['caret_command64 -file-convert -cs2fsc ' sourceedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.coord.gii source/' hemlow 'h.curv.asc']);
% system(['caret_command64 -file-convert -cs2fsc ' targetedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.coord.gii target/' hemlow 'h.curv.asc']);

system(['caret_command64 -file-convert -cs2fsc ' sourceedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.coord.gii source/' hemlow 'h.inflated.H.asc']);
system(['mris_convert -c source/' hemlow 'h.inflated.H.asc source/' hemlow 'h.sphere.asc source/' hemlow 'h.inflated.H'])
copyfile(['source/' hemlow 'h.inflated.H'],['source/' hemlow 'h.sulc'])


system(['caret_command64 -file-convert -cs2fsc ' targetedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.coord.gii target/' hemlow 'h.inflated.H.asc']);
system(['mris_convert -c target/' hemlow 'h.inflated.H.asc target/' hemlow 'h.sphere.asc target/' hemlow 'h.inflated.H'])
copyfile(['target/' hemlow 'h.inflated.H'],['target/' hemlow 'h.sulc'])

system(['caret_command64 -file-convert -sc -is GS ' whitefile ' -os FSS source/' hemlow 'h.smoothwm.asc'])
system(['mris_convert source/' hemlow 'h.smoothwm.asc source/' hemlow 'h.smoothwm'])
copyfile(['source/' hemlow 'h.smoothwm'],['target/' hemlow 'h.smoothwm'])

system(['mris_register -1 -norot -nocurv -nonorm source/' hemlow 'h.sphere.asc target/' hemlow 'h.sphere.asc source/' hemlow 'h.sphere.asc.reg']);

system(['caret_command64 -file-convert -sc -is FSS source/' hemlow 'h.sphere.asc.reg -os CARET Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii SPHERICAL CLOSED -struct ' fullhem])

%system(['caret_command64 -surface-sphere-project-unproject Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii  Edgereg.' hem '.reg_LR.coord.gii ' coordfile ' Edgereg.' hem '.def_sphere.coord.gii ' topofile]);
%system(['caret_command64 -deformation-map-create SPHERE ' coordfile ' ' topofile ' Edgereg.' hem '.reg_LR.coord.gii Edgereg.' hem '.topo.gii Edgereg.' hem '.deform_map']);
system(['caret_command64 -deformation-map-create SPHERE ' coordfile ' ' topofile ' Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii Edgereg.' hem '.deform_map']);

system(['caret_command64 -deformation-map-apply Edgereg.' hem '.deform_map METRIC_AVERAGE_TILE ' sourceedgemap ' ' sourceedgemap(1:end-9) '_deformed.func.gii'])
system(['caret_command64 -deformation-map-apply Edgereg.' hem '.deform_map METRIC_NEAREST_NODE ' parcels ' ' parcels(1:end-9) '_deformed.func.gii'])
