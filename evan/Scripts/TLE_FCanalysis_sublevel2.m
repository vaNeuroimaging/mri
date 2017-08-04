groups = {'CTL','TLE'};

mode_segmentation = read_4dfpimg_HCP('/data/cn4/evan/ROIs/mode_subcortical_label_LR_MNI_333.4dfp.img');
lhipsegmentation = (mode_segmentation==17);
rhipsegmentation = (mode_segmentation==53);

for groupnum = 1:length(groups)
    group = groups{groupnum};

tmasklist = ['/data/cn4/scratch/tunde/LUIGI/FCPROCESS_V4/' group '/NEW_TMASKLIST_mod.txt'];

[subjects tmasks] = textread(tmasklist,'%s%s');

outfolder = '/data/cn4/evan/TLE/';

load(['/data/cn4/scratch/tunde/LUIGI/EVAN/FCPROCESS_V4/' group '/QC.mat'])

etype = 'littleendian';

for s = 1:length(subjects)
    
    subname = subjects{s};
    disp(subname)
    
  
    anatavename = ['/data/cn4/scratch/tunde/LUIGI/EVAN/' group '/' subname '/atlas/' subname '_anat_ave_MNI_333.4dfp.img'];
    anataves(:,s) = read_4dfpimg_HCP(anatavename);
    %[voxelsize frames etype] = read_4dfpifh_HCP([anatavename(1:end-3) 'ifh']);
    
    dataname = ['/data/cn4/scratch/tunde/LUIGI/EVAN/FCPROCESS_V4/' group '/' subname '/' subname '_333_zmdt_resid_ntrpl_bpss_zmdt_g7.4dfp.img'];
    %[voxelsize frames etype] = read_4dfpifh([dataname(1:end-3) 'ifh']);
    data = read_4dfpimg_HCP(dataname);
    tmask = load(tmasks{s});
    data = data(:,logical(tmask));
    
    lhip_timecourse = mean(data(logical(lhipsegmentation),:),1);
    lhip_connectivity = FisherTransform(paircorr_mod(data',lhip_timecourse'));
    write_4dfpimg(lhip_connectivity,[outfolder '/' subname '_lHip.4dfp.img'],etype);
    write_4dfpifh_333_MNI([outfolder '/' subname '_lHip.4dfp.ifh'],1,etype);
    
    
    rhip_timecourse = mean(data(logical(rhipsegmentation),:),1);
    rhip_connectivity = FisherTransform(paircorr_mod(data',rhip_timecourse'));
    write_4dfpimg(rhip_connectivity,[outfolder '/' subname '_rHip.4dfp.img'],etype);
    write_4dfpifh_333_MNI([outfolder '/' subname '_rHip.4dfp.ifh'],1,etype);
    
     
%     timecourses_264 = QC(s).tc(:,:,end);
%     timecourses_264 = timecourses_264(logical(tmask),:);
%     matrix_264(:,:,s) = FisherTransform(paircorr_mod(timecourses_264));
%     matrix = matrix_264(:,:,s);
%     save([outfolder subname '_264.mat'],'matrix');
    
    

        
end
 
% write_4dfpimg(anataves,[outfolder group '_anat_aves.4dfp.img'],etype);
% write_4dfpifh_333_MNI([outfolder group '_anat_aves.4dfp.ifh'],length(subjects),etype)
% system(['niftigz_4dfp -n ' outfolder group '_anat_aves.4dfp.img ' outfolder '/' group '_anat_aves'])
% 
% write_4dfpimg(mean(anataves,2),[outfolder group '_mean_anat_aves.4dfp.img'],etype);
% write_4dfpifh_333_MNI([outfolder group '_mean_anat_aves.4dfp.ifh'],1,etype)
% system(['niftigz_4dfp -n ' outfolder group '_mean_anat_aves.4dfp.img ' outfolder '/' group '_mean_anat_aves'])
% 
% write_4dfpimg(std(anataves,[],2),[outfolder group '_std_anat_aves.4dfp.img'],etype);
% write_4dfpifh_333_MNI([outfolder group '_std_anat_aves.4dfp.ifh'],1,etype)
% system(['niftigz_4dfp -n ' outfolder group '_std_anat_aves.4dfp.img ' outfolder '/' group '_std_anat_aves'])

clear anataves QC matrix_264

end
