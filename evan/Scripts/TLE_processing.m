groups = {'CTL','TLE'};
template = '/data/cn4/evan/ROIs/EPI.nii';

for groupnum = 1:length(groups);
    
    group = groups{groupnum};
    
    tmasklist = ['/data/cn4/scratch/tunde/LUIGI/FCPROCESS_V4/' group '/NEW_TMASKLIST_mod.txt'];
    
    [subjects tmasks] = textread(tmasklist,'%s%s');
    
    for s = 1:length(subjects)
        
        outfolder = ['/data/cn4/scratch/tunde/LUIGI/EVAN/' group '/' subjects{s} '/'];
        
%         if ~exist(outfolder)
%             
             disp([group ' group, subject ' num2str(s)])
%             
%             
%             
%             mkdir(outfolder)
%             cd(outfolder)
%             
%             sub_anatave_4dfp = ['/data/cn4/scratch/tunde/LUIGI/' group '/' subjects{s} '/atlas/' subjects{s} '_anat_ave.4dfp.img'];
%             
%             system(['niftigz_4dfp -n ' sub_anatave_4dfp ' ' outfolder 'anat_ave']);
%             
%             system(['flirt -in ' outfolder 'anat_ave.nii.gz -ref ' template ' -out ' outfolder 'anat_ave_222.nii.gz -omat ' outfolder 'anat_ave2EPItemplate.mat'])
%             
%             subdata1_4dfp = ['/data/cn4/scratch/tunde/LUIGI/' group '/' subjects{s} '/bold1/' subjects{s} '_b1_xr3d.4dfp.img'];
%             subdata2_4dfp = ['/data/cn4/scratch/tunde/LUIGI/' group '/' subjects{s} '/bold2/' subjects{s} '_b2_xr3d.4dfp.img'];
%             system(['niftigz_4dfp -n ' subdata1_4dfp ' ' outfolder 'bold1']);
%             system(['niftigz_4dfp -n ' subdata2_4dfp ' ' outfolder 'bold2']);
%             
%             system(['flirt -in ' outfolder 'bold1.nii.gz -ref ' template ' -applyxfm -init ' outfolder 'anat_ave2EPItemplate.mat -out ' outfolder 'bold1_222.nii.gz'])
%             system(['flirt -in ' outfolder 'bold2.nii.gz -ref ' template ' -applyxfm -init ' outfolder 'anat_ave2EPItemplate.mat -out ' outfolder 'bold2_222.nii.gz'])
%             
%             mkdir([outfolder '/atlas'])
%             mkdir([outfolder '/bold1'])
%             mkdir([outfolder '/bold2'])
%             
%             system(['niftigz_4dfp -4 ' outfolder 'anat_ave_222.nii.gz ' outfolder 'atlas/' subjects{s} '_anat_ave_MNI_222.4dfp.img']);
             system(['t4img_4dfp none ' outfolder 'atlas/' subjects{s} '_anat_ave_MNI_222.4dfp.img  ' outfolder 'atlas/' subjects{s} '_anat_ave_MNI_333.4dfp.img -O/data/hcp-zfs/home/laumannt/HCP_Q1release/HCP_Q1_Release/seeds_7112B_to_MNI/seed_9.94_-45.52_72.63_333.4dfp.ifh'])
%             system(['niftigz_4dfp -4 ' outfolder 'bold1_222.nii.gz ' outfolder 'bold1/' subjects{s} '_b1_xr3d_222.4dfp.img']);
             system(['t4img_4dfp none ' outfolder 'bold1/' subjects{s} '_b1_xr3d_222.4dfp.img ' outfolder 'bold1/' subjects{s} '_b1_xr3d_333.4dfp.img -O/data/hcp-zfs/home/laumannt/HCP_Q1release/HCP_Q1_Release/seeds_7112B_to_MNI/seed_9.94_-45.52_72.63_333.4dfp.ifh'])
%             system(['niftigz_4dfp -4 ' outfolder 'bold2_222.nii.gz ' outfolder 'bold2/' subjects{s} '_b2_xr3d_222.4dfp.img']);
             system(['t4img_4dfp none ' outfolder 'bold2/' subjects{s} '_b2_xr3d_222.4dfp.img ' outfolder 'bold2/' subjects{s} '_b2_xr3d_333.4dfp.img -O/data/hcp-zfs/home/laumannt/HCP_Q1release/HCP_Q1_Release/seeds_7112B_to_MNI/seed_9.94_-45.52_72.63_333.4dfp.ifh'])
%             
%             delete([outfolder '*.nii.gz'])
%             
%         end
        
        
    end
end