%% Smooth volume with ROI, 120

%cohortdir = '/data/cn4/laumannt/fcMapping_redux/';
%cohs = {'C1';'C2';'C3'};
%cohs = {'AllC'};   
TR = 2.5;
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
hemnamelow = {'left';'right'};
smooth = 2.55;
processdir = '/data/cn5/selfRegulation/V4Process_nosmooth/';
outputdir = [processdir '/final_output_wROIsmooth/'];
mkdir(outputdir);
[subjects ign] = textread([processdir '/Finaltmasklist_120.txt'],'%s%s');
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
maskdir = '/data/cn5/selfRegulation/V4Process_nosmooth/subcortical_mask';
for s = 1:length(subjects)
    
    subfuncdir = [processdir '/final_output/' subjects{s}];
    subfuncvol = [subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt'];
    system(['cp ' subfuncdir '/' subfuncvol '.4dfp.* ' outputdir])
    system(['nifti_4dfp -n ' outputdir '/' subfuncvol '.4dfp.img ' outputdir '/' subfuncvol '.nii'])
    disp(['Processing subject #' num2str(s)])
       
    cd(outputdir)
    system([workbenchdir '/wb_command -volume-smoothing ' subfuncvol '.nii ' num2str(smooth) ' ' subfuncvol '_smooth2.55.nii -roi ' maskdir '/mode_subcortical_label_LR_333.nii'])

    system(['rm ' outputdir '/' subfuncvol '.4dfp.*'])
    system(['rm ' outputdir '/' subfuncvol '.nii'])
end
