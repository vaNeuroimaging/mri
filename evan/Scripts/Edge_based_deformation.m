function Edge_based_deformation(sourceedgemap,targetedgemap,parcels,hem)

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

coordfile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.coord.gii'];
topofile = ['/data/cn4/evan/fsaverage_LR32k/Conte69.' hem '.32k_fs_LR.topo.gii'];
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



$Caret5_Command -surface-sphere-project-unproject (sphere in: nativecoord nativetopo) "$AtlasSpaceFolder"/"$NativeFolder"/"$Subject"."$Hemisphere".sphere.reg.reg_LR.native.coord.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".sphere.164k_fs_"$Hemisphere".coord.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".def_sphere.164k_fs_"$Hemisphere".coord.gii "$AtlasSpaceFolder"/fsaverage/"$Subject"."$Hemisphere".164k_fs_"$Hemisphere".topo.gii


%system(['caret_command64 -surface-sphere-project-unproject Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii  Edgereg.' hem '.reg_LR.coord.gii ' coordfile ' Edgereg.' hem '.def_sphere.coord.gii ' topofile]);
%system(['caret_command64 -deformation-map-create SPHERE ' coordfile ' ' topofile ' Edgereg.' hem '.reg_LR.coord.gii Edgereg.' hem '.topo.gii Edgereg.' hem '.deform_map']);
system(['caret_command64 -deformation-map-create SPHERE ' coordfile ' ' topofile ' Edgereg.' hem '.coord.gii Edgereg.' hem '.topo.gii Edgereg.' hem '.deform_map']);

system(['caret_command64 -deformation-map-apply Edgereg.' hem '.deform_map METRIC_AVERAGE_TILE ' sourceedgemap ' ' sourceedgemap(1:end-9) '_deformed.func.gii'])
system(['caret_command64 -deformation-map-apply Edgereg.' hem '.deform_map METRIC_NEAREST_NODE ' parcels ' ' parcels(1:end-9) '_deformed.func.gii'])
