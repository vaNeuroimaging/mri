group = 'CTL';

tmasklist = ['/data/cn4/scratch/tunde/LUIGI/FCPROCESS_V4/' group '/NEW_TMASKLIST_mod.txt'];

[subjects tmasks] = textread(tmasklist,'%s%s');

outfolder = '/data/cn4/evan/TLE/on_CTL_TLE/';

% load(['/data/cn4/scratch/tunde/LUIGI/FCPROCESS_V4/' group '/final_process/QC.mat'])
% 
% 
% power_coords = load('/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/BB264_coords.txt');
% laum_coords = load('/data/cn4/laumannt/longRestingState/group_modules/264_consensus_hemsort_7112B.txt');
% clear ind
% for i = 1:264
%     
%        for t = 1:264
%                   if isequal(laum_coords(i,:),power_coords(t,:))
%                                  ind(i) = t;
%                   end
%        end
% end
clear hipdistribution
for s = 1:length(subjects)
    
    subname = subjects{s};
    disp(subname)
    
     
    %segfilefolder = ['/data/cn4/segmentation/freesurfer5_supercomputer/' subname '/mri/'];
    segfilefolder = ['/data/cn4/segmentation/freesurfer5_supercomputer/EPILEPSY_luigi_maccotta/on_CTL_TLE/' subname '/mri/'];
    system(['mri_convert -rl ' segfilefolder 'rawavg.mgz ' segfilefolder 'aparc.a2009s+aseg.mgz ' outfolder '/segmentations/Temp.nii'])
    system(['nifti_4dfp -4 ' outfolder '/segmentations/Temp.nii ' outfolder '/segmentations/' subname '.4dfp.img'])
    system(['t4img_4dfp none ' outfolder '/segmentations/' subname '.4dfp.img ' outfolder '/segmentations/' subname '_333.4dfp.img -0333'])
    
    
    
    segmentation = read_4dfpimg([outfolder '/segmentations/' subname '_333.4dfp.img']);
    [voxelsize frames etype] = read_4dfpifh([outfolder '/segmentations/' subname '_333.4dfp.ifh']);
    
    write_4dfpimg(single(segmentation==17),[outfolder '/segmentations/' subname '_lHip_333.4dfp.img'],etype);
    write_4dfpifh([outfolder '/segmentations/' subname '_lHip_333.4dfp.ifh'],1,etype);
    
    system(['cluster_4dfp ' outfolder '/segmentations/' subname '_lHip_333.4dfp.img -n100 -t.5'])
    
    lhipsegmentation = read_4dfpimg([outfolder '/segmentations/' subname '_lHip_333_clus.4dfp.img']);
    
     hipdistribution(:,s) = single(logical(lhipsegmentation));
     
    
     
     write_4dfpimg(single(segmentation==53),[outfolder '/segmentations/' subname '_rHip_333.4dfp.img'],etype);
    write_4dfpifh([outfolder '/segmentations/' subname '_rHip_333.4dfp.ifh'],1,etype);
    
    system(['cluster_4dfp ' outfolder '/segmentations/' subname '_rHip_333.4dfp.img -n100 -t.5'])
    
    rhipsegmentation = read_4dfpimg([outfolder '/segmentations/' subname '_rHip_333_clus.4dfp.img']);
    
    hipdistribution(:,s) = hipdistribution(:,s) + single(logical(rhipsegmentation));
     
    
    
    
    
    
%     dataname = ['/data/cn4/scratch/tunde/LUIGI/FCPROCESS_V4/' group '/final_process/' subname '/' subname '_333_zmdt_resid_ntrpl_bpss_zmdt_g7.4dfp.img'];
%     %[voxelsize frames etype] = read_4dfpifh([dataname(1:end-3) 'ifh']);
%     data = read_4dfpimg(dataname);
%     tmask = load(tmasks{s});
%     data = data(:,logical(tmask));
%     
%        
%     lhip_timecourse = mean(data(logical(lhipsegmentation),:),1);
%     lhip_connectivity(:,s) = FisherTransform(paircorr_mod(data',lhip_timecourse'));
%     write_4dfpimg(lhip_connectivity(:,s),[outfolder '/' subname '_lHip.4dfp.img'],etype);
%     write_4dfpifh([outfolder '/' subname '_lHip.4dfp.ifh'],frames,etype);
%     
%     
%        
%     rhip_timecourse = mean(data(logical(rhipsegmentation),:),1);
%     rhip_connectivity(:,s) = FisherTransform(paircorr_mod(data',rhip_timecourse'));
%     write_4dfpimg(rhip_connectivity(:,s),[outfolder '/' subname '_rHip.4dfp.img'],etype);
%     write_4dfpifh([outfolder '/' subname '_rHip.4dfp.ifh'],frames,etype);
     
     
     
     
     
    
    
%     
%     timecourses_264 = QC(s).tc(:,:,end);
%     timecourses_264 = timecourses_264(logical(tmask),ind);
%     matrix_264(:,:,s) = FisherTransform(paircorr_mod(timecourses_264));
%     matrix = matrix_264(:,:,s);
%     save([subname '_264.mat'],'matrix');
    
    

        
end
 
%  write_4dfpimg(mean(lhip_connectivity,2),[outfolder '/' group '_Mean_lHip.4dfp.img'],etype);
%  write_4dfpifh([outfolder '/' group '_Mean_lHip.4dfp.ifh'],frames,etype);
%  
%  write_4dfpimg(mean(rhip_connectivity,2),[outfolder '/' group '_Mean_rHip.4dfp.img'],etype);
%  write_4dfpifh([outfolder '/' group '_Mean_rHip.4dfp.ifh'],frames,etype);
%  
%  matrix = mean(matrix_264,3);
%  save([group '_mean_264.mat'],'matrix');
%  
%  figure_corrmat_powernetwork(matrix)
 
 write_4dfpimg(mean(hipdistribution,2),[outfolder '/' group '_Hip_distribution.4dfp.img'],etype);
 write_4dfpifh([outfolder '/' group '_Hip_distribution.4dfp.ifh'],1,etype);
 system(['nifti_4dfp -n ' outfolder '/' group '_Hip_distribution.4dfp.img ' outfolder '/' group '_Hip_distribution.nii'])
 
 write_4dfpimg(hipdistribution,[outfolder '/' group '_allHip.4dfp.img'],etype);
 write_4dfpifh([outfolder '/' group '_allHip.4dfp.ifh'],length(subjects),etype);
 system(['nifti_4dfp -n ' outfolder '/' group '_allHip.4dfp.img ' outfolder '/' group '_allHip.nii'])