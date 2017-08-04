function wta_grp_cifti(ciftilist,tmasklist,cortical_map,subcort_map,regionname,corr_type,sign,varargin)
%wta_grp_cifti(ciftilist,tmasklist,cortical_map,subcort_map,regionname,corr_type,sign,varargin)
%%%
%%% djg 7/1/12 Winner-take-all method for a group of subjects
%%% 
%%% This script uses a winner-take-all approach to parcellate a structure
%%% in the brain (e.g., basal ganglia) based on fc data.
%%% Output includes a correlation matrix (partial or total correlation
%%% depending on the switch 'corr_type' you specify) for each subject 
%%% (subcortical voxels X cortical region) and a 3D matrix combining across
%%% subjects (subcortical voxels X cortical region X subject), which uses
%%% the Fisher z transformed r values.
%%%
%%% This version will combine across subjects by averaging Fisher z 
%%% correlation values, and write out heat maps for each cortical region.
%%% Heat maps are made for both the partial correlations and the total
%%% correlations (since the partial r's will be low and total r's will be
%%% better for display purposes).
%%% Then the winner-take-all approach is applied to the group correlation
%%% matrix: each subcortical voxel is assigned to the cortical region it
%%% correlates with most (note: see 'sign' switch to deal with negative
%%% correlations). The winner-take-all approach is also applied to the
%%% group z score (transformed from t stats).
%%% This script will also resample the images into 111 space for
%%% presentation purposes, and mask the image with a subcortical mask
%%% (*Note: make sure it is the one you want- scroll down to bottom of the
%%% script)
%%%
%%% 8/27/12 added creation of 3D matrix (subcortical voxels X cortical regions X subjects)
%%% to be used later for group comparisons
%%% 8/30/12 added ability to pass in a regressor to regress out signal in a certain ROI from the subcortical data
%%%
%%% Input
%%% corrfile: same corrfile used in scrubbing scripts
%%% cortical_map: 4dfp image of cortical regions
%%% subcort_map: 4dfp image of subcortical structure of interest
%%% region: specify the specific subcortial region, as this will be used in the file names (e.g., bg or thal)
%%% corr_type: partial ('part') or total ('tot) correlation
%%% sign: positive ('pos') to use only positive correlations or absolute value ('abs') to include negative correlations too
%%% *varargin: (optional) 4dfp image of region whose signal you wish to regress out of the subcortical data 
%%%
%%% USAGE:  wta_grp(corrfile,cortical_map,subcort_map,region,corr_type,sign,*varargin)
%%% USAGE:  wta_grp('/path/to/corrfile','/path/to/your_cortical_map.4dfp.img','/path/to/subcortical_mask.4dfp.img','name_your_region','part' or 'tot','pos' or 'abs',*'/path/to/regress_out.4dfp.img')
%%% USAGE:  wta_grp('/path/to/corrfile','/path/to/modified_voxelwise_map.4dfp.img','/path/to/mask.4dfp.img','bg','part','pos',*'/path/to/regress_out_mask.4dfp.img')

%% Write a log file and assign a name to the cohort
cohort = 'C1';
if isempty(varargin);
    reg = 'none';
else
    reg = varargin{1,1};
end
log=fopen('wta_grp_command_log.txt','w');
fprintf(log,'wta_grp.m log\n corrfile = %s\n cortical map = %s\n subcortical map = %s\n region = %s\n correlation type = %s\n sign = %s\n regress = %s\n', ...
    ciftilist,cortical_map,subcort_map,regionname,corr_type,sign,reg);
fclose(log);

%% Read in corrfile, cortical map, and subcortical mask

%%% Read in the subjects from the corrfile
[subjectlist ciftifiles] = textread(ciftilist,'%s%s');
[subjectlist tmask_list] = textread(tmasklist,'%s%s');
subjectlist = subjectlist(1:120);


%fid = fopen(ciftilist);
%corrfile_a = textscan(fid,'%s%s%s');
%subjectlist = corrfile_a{:,3};  
%epi_list = corrfile_a{:,1};
%tmask_list = corrfile_a{:,2};

%%% Read in the cortical map
cort_datamat = ft_read_cifti_mod(cortical_map); cort_datamat = cort_datamat.data;
%[datamat frames voxelsize] = read_4dfpimg(cortical_map); 
%cort_datamat = datamat;
%cort_frames = frames;
%cort_voxelsize = voxelsize;  
cort_ind = unique(cort_datamat);
cort_ind(cort_ind == 0) = [];

%%% Read in the subcortical mask
subcort_datamat = ft_read_cifti_mod(subcort_map); subcort_datamat = subcort_datamat.data;
%[datamat frames voxelsize] = read_4dfpimg(subcort_map); 
%subcort_datamat = datamat;
%subcort_frames = frames;
%subcort_voxelsize = voxelsize;
subcort_ind = find(subcort_datamat > 0);

%%% Read in the mask of voxels to be regressed out
if isempty(varargin);
    reg_mat = 'empty';
else
    reg_mat = ft_read_cifti_mod(varargin{1,1}); reg_mat = reg_mat.data;
    %[datamat frames voxelsize] = read_4dfpimg(varargin{1,1});
    %reg_mat = datamat;
    reg_ind = find(reg_mat > 0);
end

%%% Set up the group correlation matrices to fill in
grp_corrmat = zeros(size(subcort_ind,1),size(cort_ind,1),numel(subjectlist)); %fills in with fisher z transformed r value
tot_grp_corrmat = zeros(size(subcort_ind,1),size(cort_ind,1),numel(subjectlist));
%% Loop through each subject

for i = 1:numel(subjectlist)
    ciftifile = ciftifiles{i};
    subject = subjectlist{i};
    tmask = tmask_list{i};

    %%% Read in the subject's EPI image and tmask
    cifti_datamat = ft_read_cifti_mod(ciftifile);
    if i==1
        cifti_template = cifti_datamat; cifti_template.data = zeros(size(cifti_template.data,1),1);
    end
    cifti_datamat = cifti_datamat.data;
    tmask_col = load(tmask);
    cifti_datamasked = cifti_datamat(:,logical(tmask_col)); % epi image with tmask applied
    frame_num = size(cifti_datamasked,2);

    %%% Create the time x cortical region matrix
    clear net;
    count = 1;
    for b = 1:size(cort_ind,1)
        ind = find(cort_datamat == cort_ind(b));   %index each cortical region
        timecourse = mean(cifti_datamasked(ind,:))';   %calculate mean timecourse across all voxels in a cortical region
        net(:,count) = timecourse;        %insert timecourse for each cortical region into "net" matrix     
        count = count+1;
    end
  
    %%% Create the time x voxels matrix for the subcortical region
    subcort_timecourse = cifti_datamasked(subcort_ind,:)';
    
    %%% Regress signal from a mask out of the subcortical region data
    if strcmp(reg_mat,'empty') == 1;
    else
        mask_tc = mean(cifti_datamasked(reg_ind,:))';  %create time x mask matrix
        tot_tmask = ones(length(mask_tc),1);
        [tempimg zb]=regress_nuisance(subcort_timecourse',mask_tc,tot_tmask);
        subcort_timecourse = tempimg';
    end
    
    %%% Calculate the partial correlation matrix: subcortical voxels x cortical region controlling for the other cortical regions
    %%% OR Calculate the total correlation matrix; subcortical voxels x cortical regions with r values
    if strcmp(corr_type,'part') == 1;
        corrmat = zeros(size(subcort_timecourse,2),size(net,2));  
        for t = 1:size(net,2);
            Z = net;
            Z(:,t) = [];
            part = partialcorr(subcort_timecourse,net,Z);
            corrmat(:,t) = part(:,t); 
        end
        tot_corrmat = corr(subcort_timecourse,net);     
    elseif strcmp(corr_type,'tot') == 1;
        corrmat = corr(subcort_timecourse,net);  %%% correlate
    else
        fprintf('Error: enter ''part'' or ''tot'' for corr_type\n');
    end
    
    %%% Fisher Z transform 
    fzcorrmat = FisherTransform(corrmat);  %%% Fisher z transform correlation matrix
    fzcorrmat(isnan(fzcorrmat)) = 0;
    grp_corrmat(:,:,i) = fzcorrmat;     
    tot_fzcorrmat = FisherTransform(tot_corrmat); %%% Fisher z transform total correlation matrix
    tot_grp_corrmat(:,:,i) = tot_fzcorrmat;    
    fprintf('writing fisher z correlation matrices for %s.\n',subject);
    savefzcorr = [subject '_fz_' corr_type 'corrmat.mat'];  
    save(savefzcorr,'fzcorrmat');
    savetotfzcorr = [subject '_fz_total_corrmat.mat'];  
    save(savetotfzcorr,'tot_fzcorrmat');    
end
% Write 3D correlation matrix subcortical voxels x cortical region x subject
fprintf('Writing group 3D correlation matrices.\n');
save_grpcorrmat = [cohort '_' sign '_fz_' corr_type 'corrmat.mat'];
save(save_grpcorrmat,'grp_corrmat');
save_grptotcorrmat = [cohort '_' sign '_fz_total_corrmat.mat'];
save(save_grptotcorrmat,'tot_grp_corrmat');

%% Perform t-tests across subjects to obtain the reliability of the correlation in each voxel

ttest_grp = zeros(size(grp_corrmat,1),size(grp_corrmat,2));
for l = 1:size(grp_corrmat,2);
    voxsub = grp_corrmat(:,l,:);
    voxsub1 = squeeze(voxsub);
    voxsub2 = voxsub1';
    [h,pv,ci,stats] = ttest(voxsub2);
    ttest_netmat = stats.tstat';
    ttest_grp(:,l) = ttest_netmat;  
end
% Write group t-test matrix subcortical voxels x cortical regions
numbersubs = numel(subjectlist);
fprintf('Writing group ttest matrix.\n');
save_t_grp = [num2str(numbersubs) '_tstat_' sign '_' corr_type 'corr_grp'];
save(save_t_grp,'ttest_grp');

%% Combine across subjects 
    
% Sum the correlation matrices across subjects
fz_grp_sum = zeros(size(fzcorrmat,1),size(fzcorrmat,2));
tot_fz_grp_sum = zeros(size(fzcorrmat,1),size(fzcorrmat,2));
for m = 1:numel(subjectlist)
    subj = subjectlist{m};
    fz_indivcorr = [subj '_fz_' corr_type 'corrmat.mat'];   
    load(fz_indivcorr);
	fz_grp_sum = fz_grp_sum + fzcorrmat;
    tot_fz_indivcorr = [subj '_fz_total_corrmat.mat']; 
    load(tot_fz_indivcorr);
    tot_fz_grp_sum = tot_fz_grp_sum + tot_fzcorrmat;
end
% divide by number of subjects to get the mean partial correlation
fz_grp_fx = fz_grp_sum./(numel(subjectlist));
tot_fz_grp_fx = tot_fz_grp_sum./(numel(subjectlist));
if strcmp(sign,'abs') == 1;
    abs_fz_grp_fx = abs(fz_grp_fx); %% take absolute value to obtain the strongest correlation    
    tot_abs_fz_grp_fx = abs(tot_fz_grp_fx);
elseif strcmp(sign,'pos') == 1;
    abs_fz_grp_fx = fz_grp_fx; %% do not take absolute value (even though variable called abs...)
    tot_abs_fz_grp_fx = tot_fz_grp_fx; 
else
    fprintf('Error: enter ''abs'' or ''pos'' for sign\n');
end
numbersubs = numel(subjectlist);
fprintf('Writing group average fisher z matrices.\n');
save_fzgrpfx = [num2str(numbersubs) '_fz_' sign '_' corr_type 'corr_grp_fx'];
save(save_fzgrpfx,'abs_fz_grp_fx');
save_totfzgrpfx = [num2str(numbersubs) '_fz_' sign '_total_corr_grp_fx'];
save(save_totfzgrpfx,'tot_abs_fz_grp_fx');
    
%% Winner take all approach- subcortical voxels x wta assignment in 333 space
% fisher z matrix

[maxz,maxzinds] = max(abs_fz_grp_fx,[],2);
fz_wta_ass = cort_ind(maxzinds) .* (maxz > 0);

% for j = 1:size(abs_fz_grp_fx,1)   
% 	[val position] = max(abs_fz_grp_fx(j,:)); 
%      if val <= 0.0;                           %***** Set a threshold, e.g., r > 0.1 means do not make an assignment when r < 0.1 *****
%           fz_wta_ass(j) = 0;
%      else
% 	fz_wta_ass(j) = cort_ind(position);
%     end
% end
% fz_wta_ass = fz_wta_ass';
fprintf('Writing fisher z winner-take-all matrix.\n');
savefzwta = [num2str(numbersubs) '_fz_' sign '_' corr_type 'corr_wta_ass'];
save(savefzwta,'fz_wta_ass');

% t stat matrix

[maxt,maxtinds] = max(ttest_grp,[],2);
ttest_wta_ass = cort_ind(maxtinds) .* (maxt > 0);

% for j = 1:size(ttest_grp,1)   
% 	[val position] = max(ttest_grp(j,:)); 
%      if val < 0;                   %***** Can set a threshold of t > 0 if you want
%           ttest_wta_ass(j) = 0;       
%      else
% 	ttest_wta_ass(j) = cort_ind(position);
%     end
% end
% ttest_wta_ass = ttest_wta_ass';
fprintf('Writing t stat winner-take-all matrix.\n');
savetwta = [num2str(numbersubs) '_t_' sign '_' corr_type 'corr_wta_ass'];
save(savetwta,'ttest_wta_ass');
    
%%% arrange for writing a 4dfp image and viewing, and write 4dfp img and ifh files
fprintf('Writing winner-take-all 4dfp files in 333 space.\n');

% fisher z
newname = [num2str(numbersubs) '_fz_' sign '_' corr_type 'corr_' regionname '_wta_output'];
cifti_template.data(subcort_ind) = fz_wta_ass;
ft_write_cifti_mod(newname,cifti_template);

% newimg = zeros(size(subcort_datamat,1),1);
% newimg(subcort_ind) = fz_wta_ass;
% 
% write_4dfpimg(newimg,[newname '.img'],'bigendian');
% write_4dfpifh([newname '.ifh'],1,'bigendian');
% system(['ifh2hdr ' newname '.img']);

% t stat
newnamet = [num2str(numbersubs) '_t_' sign '_' corr_type 'corr_' regionname '_wta_output'];
cifti_template.data(subcort_ind) = ttest_wta_ass;
ft_write_cifti_mod(newnamet,cifti_template);


% newimgt = zeros(size(subcort_datamat,1),1);
% newimgt(subcort_ind) = ttest_wta_ass;
% 
% write_4dfpimg(newimgt,[newnamet '.img'],'bigendian');
% write_4dfpifh([newnamet '.ifh'],1,'bigendian');
% system(['ifh2hdr ' newnamet '.img']);

%% Write heat maps for each cortical region

% fisher z
% write 4dfp images for each cortical region, resample to 111, and read back in
%networks = zeros(6443008,size(abs_fz_grp_fx,2)); %
net_name = ['allsystems_fz_' sign '_' corr_type 'corr_' regionname '_output'];
cifti_template.data = zeros(size(cifti_template.data,1),size(abs_fz_grp_fx,2));
cifti_template.data(subcort_ind,:) = abs_fz_grp_fx;
ft_write_cifti_mod(net_name,cifti_template);

net_name = ['allsystems_t_' sign '_' corr_type 'corr_' regionname '_output'];
cifti_template.data = zeros(size(cifti_template.data,1),size(ttest_grp,2));
cifti_template.data(subcort_ind,:) = ttest_grp;
ft_write_cifti_mod(net_name,cifti_template);

%% Winner take all approach in 111 space
% fisher z
% 
% for p = 1:size(abs_fz_grp_fx,2)
% 	%net_name = ['system' num2str(p) '_fz_' sign '_' corr_type 'corr_' regionname '_output'];
% 	net_col = abs_fz_grp_fx(:,p);
% 	indiv_net = zeros(size(subcort_datamat,1),1);
% 	indiv_net(subcort_ind) = net_col;
% 	write_4dfpimg(indiv_net,[net_name '.4dfp.img'],'bigendian');
% 	write_4dfpifh([net_name '.4dfp.ifh'],1,'bigendian');
% 	system(['ifh2hdr ' net_name '.4dfp.img']);
% 	system(['t4img_4dfp none ' net_name '.4dfp ' net_name '_111.4dfp -O111']);
% 	net_img = [net_name '_111.4dfp.img'];
% 	[datamat frames voxelsize] = read_4dfpimg_HCP(net_img);
% 	networks(:,p) = datamat;
% 	new_subcort_ind = find(datamat ~= 0);
% end
% % create matrix of subcortical voxels x cortical region
% new_networks = networks(new_subcort_ind,:);
% savegrpfx111 = [num2str(numbersubs) '_fz_' sign '_' corr_type 'corr_grp_fx_111'];
% save(savegrpfx111,'new_networks');
% 
% % fisher z total correlation
% % write 4dfp images for each cortical region, resample to 111, and read back in
% tot_networks = zeros(6443008,size(tot_abs_fz_grp_fx,2)); %
% for p = 1:size(tot_abs_fz_grp_fx,2)
% 	net_name = ['system' num2str(p) '_fz_' sign '_total_corr_' regionname '_output'];
% 	net_col = tot_abs_fz_grp_fx(:,p);
% 	indiv_net = zeros(size(subcort_datamat,1),1);
% 	indiv_net(subcort_ind) = net_col;
% 	write_4dfpimg(indiv_net,[net_name '.4dfp.img'],'bigendian');
% 	write_4dfpifh([net_name '.4dfp.ifh'],1,'bigendian');
% 	system(['ifh2hdr ' net_name '.4dfp.img']);
% 	system(['t4img_4dfp none ' net_name '.4dfp ' net_name '_111.4dfp -O111']);
% 	net_img = [net_name '_111.4dfp.img'];
% 	[datamat frames voxelsize] = read_4dfpimg_HCP(net_img);
% 	tot_networks(:,p) = datamat;
% 	new_subcort_indtot = find(datamat ~= 0);
% end
% % create matrix of subcortical voxels x cortical region
% tot_new_networks = tot_networks(new_subcort_indtot,:);
% savetotgrpfx111 = [num2str(numbersubs) '_fz_' sign '_total_corr_grp_fx_111'];
% save(savetotgrpfx111,'tot_new_networks');
% 
% % t stat
% % write 4dfp images for each cortical region, *t to z transform*, resample to 111, and read back in
% networksz = zeros(6443008,size(ttest_grp,2)); %
% networkst = zeros(6443008,size(ttest_grp,2)); %
% for p = 1:size(ttest_grp,2)
% 	net_namet = ['system' num2str(p) '_t_' sign '_' corr_type 'corr_' regionname '_output'];
% 	net_colt = ttest_grp(:,p);
% 	indiv_nett = zeros(size(subcort_datamat,1),1);
% 	indiv_nett(subcort_ind) = net_colt;
% 	write_4dfpimg(indiv_nett,[net_namet '.4dfp.img'],'bigendian');
% 	write_4dfpifh([net_namet '.4dfp.ifh'],1,'bigendian');
% 	system(['ifh2hdr ' net_namet '.4dfp.img']);
%     system(['t4img_4dfp none ' net_namet '.4dfp ' net_namet '_111.4dfp -O111']);
%     net_imgt = [net_namet '_111.4dfp.img'];
% 	[datamat frames voxelsize] = read_4dfpimg_HCP(net_imgt);
% 	networkst(:,p) = datamat;
% 	new_subcort_indt = find(datamat ~= 0);   
%     % transform to z scores
%     system(['t2z_4dfp ' net_namet '.img -N' num2str(numbersubs)]);    
% 	system(['t4img_4dfp none ' net_namet '_z.4dfp ' net_namet '_z_111.4dfp -O111']);
% 	net_imgz = [net_namet '_z_111.4dfp.img'];
% 	[datamat frames voxelsize] = read_4dfpimg_HCP(net_imgz);
% 	networksz(:,p) = datamat;
% 	new_subcort_indz = abs(datamat) > .00000001;
% end
% % create z score matrix of subcortical voxels x cortical region
% new_networksz = networksz(new_subcort_indz,:);
% savegrpz111 = [num2str(numbersubs) '_z_' sign '_' corr_type 'corr_grp_fx_111'];
% save(savegrpz111,'new_networksz');
% 
% %% Winner take all approach in 111 space
% % fisher z
% fprintf('Writing winner-take-all 4dfp files in 111 space.\n');
% for q = 1:size(new_networks,1)
% 	[val position] = max(new_networks(q,:)); 
%       if val <= 0.0;                           %***** Set a threshold, e.g., r > 0.1 means do not make an assignment when r < 0.1 *****
%           fz_wta_ass_111(q) = 0;        
%       else
%          fz_wta_ass_111(q) = cort_ind(position);
%      end	
% end
% fz_wta_ass_111 = fz_wta_ass_111';
% savewta111 = [num2str(numbersubs) '_fz_' sign '_' corr_type 'corr_wta_ass_111'];
% save(savewta111,'fz_wta_ass_111');
%     
% %%% arrange for writing a 4dfp image and viewing, and write 4dfp img and
% %%% ifh files
% newimg111 = zeros(6443008,1);
% newimg111(new_subcort_ind) = fz_wta_ass_111;
% newname111 = [num2str(numbersubs) '_fz_' sign '_' corr_type 'corr_' regionname '_wta_output_111'];
% write_4dfpimg(newimg111,[newname111 '.4dfp.img'],'bigendian');
% write_4dfpifh_new([newname111 '.4dfp.ifh'],1,'bigendian',176,208,176,1,1,1);
% system(['ifh2hdr ' newname111 '.4dfp.img']);
% 
% % z score
% for qq = 1:size(new_networksz,1)
% 	[val position] = max(new_networksz(qq,:)); 
%     if val <= 0.000001;
%         zscore_wta_ass_111(qq) = 0;
%     else
%         zscore_wta_ass_111(qq) = cort_ind(position);
%     end
% end
% zscore_wta_ass_111 = zscore_wta_ass_111';
% savewta111z = [num2str(numbersubs) '_z_' sign '_' corr_type 'corr_wta_ass_111'];
% save(savewta111z,'zscore_wta_ass_111');
% %     
% % %%% arrange for writing a 4dfp image and viewing, and write 4dfp img and
% % %%% ifh files
% newimg111z = zeros(6443008,1);
% newimg111z(new_subcort_indz) = zscore_wta_ass_111;
% newname111z = [num2str(numbersubs) '_z_' sign '_' corr_type 'corr_' regionname '_wta_output_111.4dfp'];
% write_4dfpimg(newimg111z,[newname111z '.img'],'bigendian');
% write_4dfpifh_new([newname111z '.ifh'],1,'bigendian',176,208,176,1,1,1);
% system(['ifh2hdr ' newname111z '.img']);
% 
% %% Mask the winner-take-all map by the subcortical structure
% bg_mask = ['/data/cn4/dgreene/BG_WTA/masks/180_' regionname '_intersect135_binarymask.4dfp.img'];
% %bg_mask = ['/data/cn4/dgreene/BG_WTA/masks/180_' region '_intersect45_binarymask.4dfp.img'];
% system(['maskimg_4dfp ' newname111 '.4dfp.img ' bg_mask ' ' newname111 '_mask']);
% system(['maskimg_4dfp ' newname111z '.4dfp.img ' bg_mask ' ' newname111z '_mask']);

%% Done!
fprintf('DONE!\n')
end