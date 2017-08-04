function fsl_preprocess_v2(folder,BOLDrunbasenames,MPRAGEbase,T2base,TR)

template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm.nii.gz';


%find T1 and T2 files
MPRAGEfile = dir([folder '/' MPRAGEbase '*.nii.gz']);
MPRAGEfile = [folder '/' MPRAGEfile(1).name(1:end-7)];
movefile(MPRAGEfile,[folder '/T1.nii.gz'])
MPRAGEfile = [folder '/T1.nii.gz'];

T2file = dir([folder '/' T2base '*.nii.gz']);
T2file = [folder '/' T2file(1).name(1:end-7)];
movefile(T2file, [folder '/T2.nii.gz'])
T2file = [folder '/T2.nii.gz'];

%reorient, crop, and bias-correct T1 and T2 images
system(['fsl_anat -i ' MPRAGEfile ' -o struct --clobber --noreg --nononlinreg --noseg --nosubcortseg'])
copyfile([folder '/struct.anat/T1_biascorr.nii.gz'],[folder '/T1_biascorr.nii.gz'])
MPRAGEfile = [folder '/T1_biascorr.nii.gz'];

system(['fsl_anat -i ' T2file ' -o struct -t T2 --clobber --noreg --nononlinreg --noseg --nosubcortseg'])
copyfile([folder '/struct.anat/T2_biascorr.nii.gz'],[folder '/T2_biascorr.nii.gz'])
T2file = [folder '/T2_biascorr.nii.gz'];

cd(folder)

%brain-extract T1 and T2 images
system(['bet2 ' MPRAGEfile ' ' MPRAGEfile '_bet']);
system(['bet2 ' T2file ' ' T2file '_bet']);

%register T1 to MNI, T2 to T1, and concatenate
system(['flirt -in ' MPRAGEfile '_bet -ref ' template ' -omat T1_2MNI.mat']);
system(['flirt -in ' MPRAGEfile '_bet -ref ' template ' -applyxfm -init T1_2MNI.mat -out ' MPRAGEfile '_bet_MNI']);
system(['flirt -in ' T2file '_bet -ref ' MPRAGEfile '_bet -omat T2_2T1.mat']);
system(['concat_xfm -omat T2_2MNI.mat -concat T1_2MNI.mat T2_2T1.mat'])


for runnum = 1:length(BOLDrunbasenames)
    BOLDfile = dir([folder '/' BOLDrunbasenames{runnum} '*.nii.gz']);
    BOLDfile = [folder '/' BOLDfile.name(1:end-7)];
    
    %motion-correct BOLD images
    system(['mcflirt -in ' BOLDfile ' -plots ']);
    
    %slice-time correct BOLD images
    system(['slicetimer -i ' BOLDfile '_mcf -o ' BOLDfile '_mcf_st -r ' num2str(TR) ' --odd'])
    
    %register BOLD to T2 and concatenate through to MNI
    system(['flirt -in ' BOLDfile '_mcf_st -ref '  T2file '_bet -omat ' BOLDfile '_2T2.mat']);
    system(['concat_xfm -omat ' BOLDfile '_2MNI.mat -concat T2_2MNI.mat ' BOLDfile '_2T2.mat'])
    
    %apply concatenated registration
    system(['flirt -in ' BOLDfile '_mcf -ref ' template ' -applyisoxfm 3 -init ' BOLDfile '_2MNI.mat -out ' BOLDfile '_mcf_st_MNI']);
    
end
    
    
    