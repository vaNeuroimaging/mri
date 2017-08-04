function Edge_based_deformation_SD(sourceedgemap,targetedgemap,parcels,hem)
%Edge_based_deformation_SD(sourceedgemap,targetedgemap,[parcels],hem)



thisfolder = pwd;

cd /data/cn4/evan/Scripts/SDv1.5.1-svn593/
add_all_paths
cd(thisfolder)

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

mkdir([outfolder '/source/surf/'])
mkdir([outfolder '/target/surf/'])

coordfile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.coord.gii'];
topofile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.32k_fs_LR.topo.gii'];
sulcfile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.sulc.32k_fs_LR_' hem '.func.gii'];
whitefile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.white.32k_fs_LR.surf.gii'];
spherefile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii'];

system(['caret_command64 -file-convert -sc -is CARET ' coordfile ' ' topofile ' -os FSS ' outfolder '/source/surf/' hemlow 'h.sphere.asc'])
system(['mris_convert ' outfolder '/source/surf/' hemlow 'h.sphere.asc ' outfolder '/source/surf/' hemlow 'h.sphere'])
copyfile(['' outfolder '/source/surf/' hemlow 'h.sphere'],['' outfolder '/target/surf/' hemlow 'h.sphere'])

% system(['caret_command64 -file-convert -cs2fsc ' sulcfile ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.coord.gii ' outfolder '/source/surf/' hemlow 'h.sulc.asc']);
% system(['mris_convert -c ' outfolder '/source/surf/' hemlow 'h.sulc.asc ' outfolder '/source/surf/' hemlow 'h.sphere.asc ' outfolder '/source/surf/' hemlow 'h.sulc'])
% copyfile(['' outfolder '/source/surf/' hemlow 'h.sulc'],['' outfolder '/target/surf/' hemlow 'h.sulc'])


% system(['caret_command64 -file-convert -cs2fsc ' sourceedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.coord.gii ' outfolder '/source/surf/' hemlow 'h.curv.asc']);
% system(['caret_command64 -file-convert -cs2fsc ' targetedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.coord.gii ' outfolder '/target/surf/' hemlow 'h.curv.asc']);

system(['caret_command64 -file-convert -cs2fsc ' sourceedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.coord.gii ' outfolder '/source/surf/' hemlow 'h.inflated.H.asc']);
system(['mris_convert -c ' outfolder '/source/surf/' hemlow 'h.inflated.H.asc ' outfolder '/source/surf/' hemlow 'h.sphere ' outfolder '/source/surf/' hemlow 'h.inflated.H'])
copyfile(['' outfolder '/source/surf/' hemlow 'h.inflated.H'],['' outfolder '/source/surf/' hemlow 'h.sulc'])
copyfile(['' outfolder '/source/surf/' hemlow 'h.inflated.H'],['' outfolder '/source/surf/' hemlow 'h.curv'])


system(['caret_command64 -file-convert -cs2fsc ' targetedgemap ' 1 /data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.coord.gii ' outfolder '/target/surf/' hemlow 'h.inflated.H.asc']);
system(['mris_convert -c ' outfolder '/target/surf/' hemlow 'h.inflated.H.asc ' outfolder '/target/surf/' hemlow 'h.sphere ' outfolder '/target/surf/' hemlow 'h.inflated.H'])
copyfile(['' outfolder '/target/surf/' hemlow 'h.inflated.H'],['' outfolder '/target/surf/' hemlow 'h.sulc'])
copyfile(['' outfolder '/target/surf/' hemlow 'h.inflated.H'],['' outfolder '/target/surf/' hemlow 'h.curv'])

system(['caret_command64 -file-convert -sc -is GS ' whitefile ' -os FSS ' outfolder '/source/surf/' hemlow 'h.smoothwm.asc'])
system(['mris_convert ' outfolder '/source/surf/' hemlow 'h.smoothwm.asc ' outfolder '/source/surf/' hemlow 'h.smoothwm'])
copyfile(['' outfolder '/source/surf/' hemlow 'h.smoothwm'],['' outfolder '/target/surf/' hemlow 'h.smoothwm'])

system(['mris_register -1 -norot -nocurv -nonorm ' outfolder '/source/surf/' hemlow 'h.sphere.asc ' outfolder '/target/surf/' hemlow 'h.sphere.asc ' outfolder '/source/surf/' hemlow 'h.sphere.asc.reg']);
cd /data/cn4/evan/Scripts/SDv1.5.1-svn593/SphericalDemons/freesurfer/
mris_SD_pairwise_register(['' outfolder '/source/surf/' hemlow 'h.sphere'],['' outfolder '/target/surf/' hemlow 'h.sphere'],['' outfolder '/source/surf/' hemlow 'h.sphere.asc.reg']);
cd(outfolder)

system(['caret_command64 -file-convert -sc -is FSS ' outfolder '/source/surf/' hemlow 'h.sphere.asc.reg -os CARET Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii SPHERICAL CLOSED -struct ' fullhem])
system(['caret_command64 -file-convert -sc -is CARET Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii -os GS Edgereg.' hem '.surf.gii'])
%system(['caret_command64 -file-convert -sc -is CARET Edgereg.' hem '.coord.gii ' topofile ' -os GS Edgereg.' hem '.def_surf.gii'])
system(['caret_command64 -file-convert -sc -is CARET ' coordfile ' Edgereg.' hem '.topo.gii -os GS Edgereg.' hem '.def_surf.gii'])

% 
% %$Caret5_Command -surface-sphere-project-unproject (sphere in: nativecoord nativetopo) "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.coord.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".sphere.164k_fs_"$Hemisphere".coord.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere.164k_fs_"$Hemisphere".coord.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".164k_fs_"$Hemisphere".topo.gii
% 
% 

%system(['/data/cn4/evan/workbench/bin_rh_linux64/wb_command -surface-sphere-project-unproject Edgereg.' hem '.surf.gii ' spherefile ' Edgereg.' hem '.def_surf.gii Edgereg.' hem '.reg_LR.surf.gii'])
system(['/data/cn4/evan/workbench/bin_linux64/wb_command -surface-sphere-project-unproject ' spherefile ' Edgereg.' hem '.surf.gii Edgereg.' hem '.def_surf.gii Edgereg.' hem '.reg_LR.surf.gii'])
system(['caret_command64 -file-convert -sc -is GS Edgereg.' hem '.reg_LR.surf.gii -os CARET Edgereg.' hem '.reg_LR.coord.gii Edgereg.' hem '.reg_LR.topo.gii SPHERICAL CLOSED -struct ' fullhem])

% %system(['caret_command64 -surface-sphere-project-unproject Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii  Edgereg.' hem '.reg_LR.coord.gii ' coordfile ' Edgereg.' hem '.def_sphere.coord.gii ' topofile]);

system(['caret_command64 -deformation-map-create SPHERE ' coordfile ' ' topofile ' Edgereg.' hem '.reg_LR.coord.gii Edgereg.' hem '.topo.gii Edgereg.' hem '.deform_map']);
% system(['caret_command64 -deformation-map-create SPHERE ' coordfile ' ' topofile ' Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii Edgereg.' hem '.deform_map']);
% 
system(['caret_command64 -deformation-map-apply Edgereg.' hem '.deform_map METRIC_AVERAGE_TILE ' sourceedgemap ' ' sourceedgemap(1:end-9) '_deformed.func.gii'])
if ~isempty(parcels)
    system(['caret_command64 -deformation-map-apply Edgereg.' hem '.deform_map METRIC_NEAREST_NODE ' parcels ' ' parcels(1:end-9) '_deformed.func.gii'])
end

delete('Edgereg*')
delete('spec_file_128.spec')
delete('error.log')
rmdir('source')
rmdir('target')
