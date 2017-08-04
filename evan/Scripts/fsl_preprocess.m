function fsl_preprocess(folder,BOLDrunbasenames,MPRAGEbase)

template = '/usr/local/fsl/data/standard/MNI152lin_T1_1mm.nii.gz';

MPRAGEfile = dir([folder '/' MPRAGEbase '*.nii.gz']);
MPRAGEfile = [folder '/' MPRAGEfile(1).name(1:end-7)];
movefile(MPRAGEfile,[folder '/T1.nii.gz'])
MPRAGEfile = [folder '/T1.nii.gz'];

system(['fsl_anat -i ' MPRAGEfile ' --clobber --noreg ----nononlinreg --noseg --nosubcortseg'])
copyfile([folder '/T1.anat/T1_biascorr.nii.gz'],[folder '/T1_biascorr.nii.gz'])
MPRAGEfile = [folder '/T1_biascorr.nii.gz'];

T2file = dir([folder '/' T2base '*.nii.gz']);
T2file = [folder '/' T2file(1).name(1:end-7)];
movefile(T2file, [folder '/T2.nii.gz'])
T2file = [folder '/T2.nii.gz'];

cd(folder)

system(['bet2 ' MPRAGEfile ' ' MPRAGEfile '_bet']);
system(['flirt -in ' MPRAGEfile '_bet -ref ' template ' -omat T12MNI.mat']);
system(['flirt -in ' MPRAGEfile '_bet -ref ' template ' -applyxfm -init T12MNI.mat -out ' MPRAGEfile '_bet_MNI']);
 
system(['bet2 ' T2file ' ' T2file '_bet']);
system(['flirt -in ' T2file '_bet -ref ' MPRAGEfile '_bet -omat T2_2T1.mat']);

%T12MNI = load('T12MNI.mat','-ascii');
%T2_2T1 = load('T2_2T1.mat','-ascii');

for runnum = 1:length(BOLDrunbasenames)
    BOLDfile = dir([folder '/' BOLDrunbasenames{runnum} '*.nii.gz']);
    BOLDfile = [folder '/' BOLDfile.name(1:end-7)];
    
    system(['mcflirt -in ' BOLDfile ' -plots ']);
    
    %system(['flirt -in ' BOLDfile '_mcf -ref '  T2file '_bet -omat ' BOLDfile '_2T2.mat']);
    
    %BOLD2T2 = load([BOLDfile '_2T2.mat'],'-ascii');
    %BOLD2MNI = BOLD2T2 * T2_2T1 * T12MNI;
    %save([BOLDfile '_2MNI.mat'],'BOLD2MNI','-ascii')
    
    %system(['flirt -in ' BOLDfile '_mcf -ref ' template ' -applyisoxfm 3 -init ' BOLDfile '_2MNI.mat -out ' BOLDfile '_mcf_MNI']);
    
    
    system(['epi_reg --epi=' BOLDfile '_mcf --t1=' MPRAGEfile ' --t1brain=' MPRAGEfile '_bet --out=' BOLDfile '_mcf_T1']);
    system(['flirt -in ' BOLDfile '_mcf_T1 -ref ' template ' -applyisoxfm 3 -init T12MNI.mat -out ' BOLDfile '_mcf_T1_MNI']);
    
    
end
    
    
    