function systemresult = PostFreeSurferPipeline_fsavg2fslr_long_singlesub(subject,T1name,FreesurferImportLocation,StudyFolder,preprocessedfolder,systemresult)
% CG: adapted .bat function with the same name to work in matlab for easier


% Paths
InputAtlasName = 'NativeVol'; %e.g., 7112b, DLBS268, MNI152 (will be used to name folders in subject's folder)
ResampleAtlasName = 'MNI';

if ~exist('systemresult')
    systemresult = cell(0,2);
end

%FreesurferImportLocation='/home/data/cn4/segmentation/freesurfer5_supercomputer/adults/'; %Input freesurfer path
%FreesurferImportLocation=['/home/data/subjects/' subject '/freesurfer/'];
%MPRImportLocation='/home/data/cn4/segmentation/freesurfer5_supercomputer/adults/'; %Input MPR path
MPRImportLocation=FreesurferImportLocation;
%StudyFolder='/home/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR'; %Path to overall study folder
%StudyFolder=['/home/data/subjects/' subject '/fs_LR/'];
mkdir(StudyFolder)
cd(StudyFolder)
%T1name = 'T1_biascorr.nii.gz';

FinalTemplateSpace=[FreesurferImportLocation '/' T1name];

DownSampleI='32000'; %Downsampled mesh (same as HCP)
DownSampleNameI='32';

%Script and Program Locations:
CaretAtlasFolder='/home/data/scripts/Cifti_creation/gagan_tools/standard_mesh_atlases'; %Copied from Git Repo
PipelineScripts='/home/data/scripts/Cifti_creation/gagan_tools/scripts'; %Copied from Git Repo
%PipelineBinaries='/home/data/scripts/Cifti_creation/gagan_tools/global/binaries'; %Copied from Git Repo
GlobalScripts='/home/data/scripts/Cifti_creation/gagan_tools/global/scripts'; %Copied from Git Repo
Caret5_Command='/usr/local/caret/bin_linux64/caret_command'; %Location of Caret5 caret_command
Caret7_Command='/usr/local/workbench/bin_rh_linux64/wb_command';%Location of Caret7 wb_command

%Image locations and names:
T1wFolder=[StudyFolder '/' InputAtlasName]; %Could try replacing "$InputAtlasName" everywhere with String of your choice, e.g. 7112bLinear
AtlasSpaceFolder=[StudyFolder '/' InputAtlasName];
NativeFolder='Native';
FreeSurferFolder=FreesurferImportLocation;
FreeSurferInput=T1name;
T1wRestoreImage=T1name;
T2wRestoreImage = T1name;
AtlasTransform=[StudyFolder '/' InputAtlasName '/zero'];
InverseAtlasTransform=[StudyFolder '/' InputAtlasName '/zero'];
AtlasSpaceT1wImage = T1name;
AtlasSpaceT2wImage = T1name;
T1wImageBrainMask='brainmask_fs'; %Name of FreeSurfer-based brain mask -- I think this gets created? GW

%Making directories and copying over relevant data (freesurfer output and mpr):
[systemresult{end+1,1},systemresult{end+1,2}] = system(['mkdir -p ' StudyFolder '/' InputAtlasName]);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['cp -R ' FreesurferImportLocation ' ' StudyFolder '/' InputAtlasName]);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['cp ' MPRImportLocation '/' T1name '.nii.gz ' StudyFolder '/' InputAtlasName '/' T1name '.nii.gz']);

%I think this stuff below is making the 'fake warpfield that is identity above? GW
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' StudyFolder '/' InputAtlasName '/' T1name ' -sub ' StudyFolder '/' InputAtlasName '/' T1name ' ' StudyFolder '/' InputAtlasName '/zero.nii.gz']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmerge -t ' StudyFolder '/' InputAtlasName '/zero_.nii.gz ' StudyFolder '/' InputAtlasName '/zero.nii.gz ' StudyFolder '/' InputAtlasName '/zero.nii.gz ' StudyFolder '/' InputAtlasName '/zero.nii.gz']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['mv -f ' StudyFolder '/' InputAtlasName '/zero_.nii.gz ' StudyFolder '/' InputAtlasName '/zero.nii.gz']);

% Run it
[systemresult{end+1,1},systemresult{end+1,2}] = system(['/home/data/scripts/Cifti_creation/FreeSurfer2CaretConvertAndRegisterNonlinear.sh ' StudyFolder ' ' subject ' ' T1wFolder ' ' AtlasSpaceFolder ' ' NativeFolder ' ' FreeSurferFolder ' ' FreeSurferInput ' ' FinalTemplateSpace ' ' T1wRestoreImage ' ' T2wRestoreImage ' ' CaretAtlasFolder ' ' DownSampleI ' ' DownSampleNameI ' ' Caret5_Command ' ' Caret7_Command ' ' AtlasTransform ' ' InverseAtlasTransform ' ' AtlasSpaceT1wImage ' ' AtlasSpaceT2wImage ' ' T1wImageBrainMask ' ' PipelineScripts ' ' GlobalScripts]);



%----------------------------------------------------------------------

ResampleFolder = [StudyFolder '/' ResampleAtlasName];
mkdir(ResampleFolder);
Atlasvol = [preprocessedfolder '/' T1name '_' ResampleAtlasName '.nii.gz'];
copyfile(Atlasvol,ResampleFolder);

[systemresult{end+1,1},systemresult{end+1,2}] = system(['cp -r ' AtlasSpaceFolder '/* ' ResampleFolder '/']);

transform_matrix = [preprocessedfolder '/T1_2' ResampleAtlasName '.mat'];
world_matrix = [preprocessedfolder '/T1_2' ResampleAtlasName '.world'];
[systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -convert-affine -from-flirt ' transform_matrix ' ' FinalTemplateSpace '.nii.gz ' Atlasvol ' -to-world ' world_matrix]);

surfaces = {'midthickness', 'white', 'pial', 'inflated', 'very_inflated'};

spaces = {'./','Native','fsaverage_LR32k'};
extensions_per_space = {'164k_fs_LR','native','32k_fs_LR'};

for i = 1:length(spaces)
    thisspace = spaces{i};
    thisextension = extensions_per_space{i};
    thisResampleFolder = [ResampleFolder '/' thisspace];
    cd(thisResampleFolder);
    [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -add-to-spec-file ' thisResampleFolder '/' subject '.' thisextension '.wb.spec INVALID ' ResampleFolder '/' T1name '_' ResampleAtlasName '.nii.gz']);
    
    %         [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' Subject '.L.sphere.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' Subject '.L.sphere.' thisextension '.surf.gii']);
    %         [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret5_command ' -surface-apply-transformation-matrix ' thisResampleFolder '/' Subject '.L.sphere.' thisextension '.coord.gii ' thisResampleFolder '/' Subject '.L.' thisextension '.topo.gii ' thisResampleFolder '/' Subject '.L.sphere.' thisextension '.coord.gii -matrix ' transform_matrix]);
    %         [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' Subject '.R.sphere.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' Subject '.R.sphere.' thisextension '.surf.gii']);
    %         [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret5_command ' -surface-apply-transformation-matrix ' thisResampleFolder '/' Subject '.R.sphere.' thisextension '.coord.gii ' thisResampleFolder '/' Subject '.R.' thisextension '.topo.gii ' thisResampleFolder '/' Subject '.R.sphere.' thisextension '.coord.gii -matrix ' transform_matrix]);
    %
    %         if strcmp(thisspace,'Native')
    %             [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' Subject '.L.sphere.reg.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' Subject '.L.sphere.reg.' thisextension '.surf.gii']);
    %             [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' Subject '.L.sphere.reg.reg_LR.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' Subject '.L.sphere.reg.reg_LR.' thisextension '.surf.gii']);
    %             [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' Subject '.R.sphere.reg.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' Subject '.R.sphere.reg.' thisextension '.surf.gii']);
    %             [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' Subject '.R.sphere.reg.reg_LR.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' Subject '.R.sphere.reg.reg_LR.' thisextension '.surf.gii']);
    %         end
    
    for j = 1:length(surfaces)
        surface = surfaces{j};
        
        [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' subject '.L.' surface '.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' subject '.L.' surface '.' thisextension '.surf.gii']);
        [systemresult{end+1,1},systemresult{end+1,2}] = system(['set matrix = `cat ' transform_matrix '`; ' Caret5_Command ' -surface-apply-transformation-matrix ' thisResampleFolder '/' subject '.L.' surface '.' thisextension '.coord.gii ' thisResampleFolder '/' subject '.L.' thisextension '.topo.gii ' thisResampleFolder '/' subject '.L.' surface '.' thisextension '.coord.gii -matrix $matrix']);
        
        [systemresult{end+1,1},systemresult{end+1,2}] = system([Caret7_Command ' -surface-apply-affine ' thisResampleFolder '/' subject '.R.' surface '.' thisextension '.surf.gii ' world_matrix ' ' thisResampleFolder '/' subject '.R.' surface '.' thisextension '.surf.gii']);
        [systemresult{end+1,1},systemresult{end+1,2}] = system(['set matrix = `cat ' transform_matrix '`; ' Caret5_Command ' -surface-apply-transformation-matrix ' thisResampleFolder '/' subject '.R.' surface '.' thisextension '.coord.gii ' thisResampleFolder '/' subject '.R.' thisextension '.topo.gii ' thisResampleFolder '/' subject '.R.' surface '.' thisextension '.coord.gii -matrix $matrix']);
        
    end
    
end

[systemresult{end+1,1},systemresult{end+1,2}] = system(['rm ' ResampleFolder '/fsaverage/*coord*']);









