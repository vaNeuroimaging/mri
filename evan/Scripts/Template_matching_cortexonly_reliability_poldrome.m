

xdistance = 30;

minclustersizemm2 = 30;

kden_thresh = .05;

ncortverts = 59412;

%-----------------------------


%load('/data/cn4/evan/RestingState/Ind_variability/Templates_Yeo.mat')
load('/data/cn4/evan/RestingState/Ind_variability/Templates_consensus.mat')
% templates{1}(:,11:12) = [];
% IDs{1}(11:12) = [];
templates{1}(:,end-1:end) = [];
IDs{1}(end-1:end) = [];

templates = templates{1};
IDs = IDs{1};

%%load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat
% distances = smartload('/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/normalwall_distmat_333/distmat_surf_geodesic_vol_euc.mat');
% distances = distances(1:ncortverts,1:ncortverts);
% 
% threshdistance = distances > xdistance;

clear distances


groupavg = ft_read_cifti_mod('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
groupavg = groupavg.data(1:ncortverts,1);

%%

values_sorted = sort(templates(:),'descend');
threshval = values_sorted(round(numel(templates) .* kden_thresh));
ThreshTemplates = templates >= threshval;
clear values_sorted

threshdistance = threshdistance(1:ncortverts,1:ncortverts);
ThreshTemplates = ThreshTemplates(1:ncortverts,:);




funcdir = '/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/FCPROCESS_SCRUBBED_meanfield_LSinterp_333_sub018_reg_FD025_mirpad/';
[subjects tmasks] = textread('/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/poldrome_allses_selected_final_TMASKLIST.txt','%s%s');


iter = 10;
testsessionvec = [1:3:42];
startpos = 51;
% Load cifti timeseries for all data
% for s = 1:length(subjects)
%     string{s} = ['Loading group, subject #: ' num2str(s)];
%     ciftiname = [funcdir '/cifti_timeseries_normalwall/' subjects{s} '_333_zmdt_resid_ntrpl_bpss_zmdt_wROI255_LR_surf_subcort.dtseries.nii'];
%     if s==1; neighbors = cifti_neighbors(ciftiname); neighbors = neighbors(1:ncortverts,:); fprintf('%s',string{s}); else fprintf([repmat('\b',1,length(string{s-1})) '%s'],string{s}); end
%     tmask = load(tmasks{s}); tmask(1:(startpos-1)) = 0;
%     temp = ft_read_cifti_mod(ciftiname);
%     ciftifiles{s} = single(temp.data(1:ncortverts,logical(tmask)));
%     
% end
% disp(' ')


overlap = zeros(iter,length(testsessionvec));


for n = 1:iter
    disp(['iter #' num2str(n)])
    
    order = randperm(length(subjects));
%     otherhalf_inds = order((length(subjects)/2):end);
%     otherhalf = [ciftifiles{otherhalf_inds}];
%     
%     %otherhalf = time_concat((this_iterind > (length(subjects)/2)),:);
%     
%     
%     %-----------------------------
%     %Calculate and threshold correlations
%     
%     
%     correlmaps = paircorr_mod(otherhalf');
%     clear otherhalf
%     correlmaps(isnan(correlmaps)) = 0;
%     
%     temp=correlmaps(triu(true(ncortverts),1));
%     v=sort(temp,'descend');
%     clear temp
%     r_sub_thresh = v(round(kden_thresh * numel(v)));
%     clear v
%     correlmaps_thresh = correlmaps > r_sub_thresh;
%     
%     clear correlmaps
%     
%     prevstring = [];
%     
%     
%     %-----------------------------
%     %Match to templates
%     
%     for templatenum = 1:size(ThreshTemplates,2);
%         string = ['Half of data: comparing to network ' num2str(IDs(templatenum))];
%         fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%         prevstring = string;
%         template_big = repmat(ThreshTemplates(:,templatenum)',ncortverts,1);
%         dice_coeffs(:,templatenum) = sum((template_big & correlmaps_thresh & threshdistance),2) ./ sum(((template_big | correlmaps_thresh) & threshdistance),2);
%     end
%     clear template_big
%     
%     [sorteddice, sorti] = sort(dice_coeffs,2,'descend');
%     sorti = sorti(:,1)';
%     dice_maxes = sorteddice(:,1);
%     clear sorteddice
%     
%     otherhalf_maps = IDs(sorti)';
%     
%     otherhalf_maps(dice_maxes==0) = 0;
%     
%     clear dice_coeffs
    
    for t = 1:length(testsessionvec)
        
        %testdata = time_concat((this_iterind <= testsessionvec(t)),:);
        testdata = [ciftifiles{order(1:t)}];
        
        correlmaps = paircorr_mod(testdata');
        clear testdata
        correlmaps(isnan(correlmaps)) = 0;
        
        temp=correlmaps(triu(true(ncortverts),1));
        v=sort(temp,'descend');
        clear temp
        r_sub_thresh = v(round(kden_thresh * numel(v)));
        clear v
        correlmaps_thresh = correlmaps > r_sub_thresh;
        
        clear correlmaps
        
        
        prevstring = [];
        
        %-----------------------------
        %Match to templates
        
        for templatenum = 1:size(ThreshTemplates,2);
            string = ['Test data from ' num2str(testsessionvec(t)) ' sessions: comparing to network ' num2str(IDs(templatenum))];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            template_big = repmat(ThreshTemplates(:,templatenum)',ncortverts,1);
            dice_coeffs(:,templatenum) = sum((template_big & correlmaps_thresh & threshdistance),2) ./ sum(((template_big | correlmaps_thresh) & threshdistance),2);
        end
        clear template_big
        
        [sorteddice, sorti] = sort(dice_coeffs,2,'descend');
        sorti = sorti(:,1)';
        dice_maxes = sorteddice(:,1);
        clear sorteddice
        
        testdata_maps = IDs(sorti)';
        
        testdata_maps(dice_maxes==0) = 0;
        
        %nonzeroinds = (testdata_maps>0) | (otherhalf_maps>0);
        
        %overlap(n,t) = nnz(testdata_maps(nonzeroinds)==otherhalf_maps(nonzeroinds)) / nnz(nonzeroinds);
        
        nonzeroinds = (testdata_maps>0) | (groupavg>0);
        
        overlap(n,t) = nnz(testdata_maps(nonzeroinds)==groupavg(nonzeroinds)) / nnz(nonzeroinds);
        
        
    end
    disp(' ')
end









%
%         %-----------------------------
%         %Size threshold
%
%
%
%         string = ['Subject ' num2str(subnums(s)) ': eliminating patches smaller than ' num2str(minclustersizemm2) 'mm2'];
%         fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%         prevstring = string;
%
%
%         hems = {'L','R'};
%         for hemnum = 1:length(hems)
%             surfaceareafiles{hemnum} = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
%             %surfaceareafiles{hemnum} = ['/data/hcp-zfs/shared-nil/adeyemob/HCP_MAIN/OTHER/fsaverage_LR32K/' subject '/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
%             if ~exist(surfaceareafiles{hemnum})
%                 surffile = [surfaceareafiles{hemnum}(1:end-22) '.surf.gii'];
%                 system(['wb_command -surface-vertex-areas ' surffile ' ' surfaceareafiles{hemnum}])
%             end
%             surfaceareas{hemnum} = gifti(surfaceareafiles{hemnum});
%             surfaceareas{hemnum} = surfaceareas{hemnum}.cdata;
%             surfaceareas{hemnum}(out_template.brainstructure((1:length(surfaceareas{hemnum})) + (length(surfaceareas{hemnum}) * (hemnum-1)))==-1) = [];
%         end
%
%         surfacearea_voxvol = [surfaceareas{1} ; surfaceareas{2} ; (ones(size(correlmaps_thresh,1) - ncortverts , 1) * voxvol)];
%
%         temp = zeros(ncortverts,1);
%         for ID = IDs(:)'
%             clustereddata = cifti_cluster_surfacearea_volume(thissub_networkconnections,ID-.5,ID+.5,minclustersizemm2,minclustersizemm2,ncortverts,surfacearea_voxvol,neighbors);
%             clustereddata(logical(clustereddata)) = ID;
%             temp = temp + clustereddata;
%         end
%
%         thissub_networkconnections = temp;
%
%
%
%         %-----------------------------
%         %Fill in removed spots
%
%
%
%         blankinds = find(thissub_networkconnections==0);
%         while ~isempty(blankinds)
%             temp = thissub_networkconnections;
%             for ind = blankinds(:)'
%                 indneighs = neighbors(ind,2:end); indneighs(isnan(indneighs)) = [];
%                 neighvals = thissub_networkconnections(indneighs); neighvals(neighvals==0) = [];
%                 if ~isempty(neighvals)
%                     temp(ind) = mode(neighvals);
%                 end
%             end
%             thissub_networkconnections = temp;
%             blankinds = find(thissub_networkconnections==0);
%         end
%
%         thissub_networkconnections(dice_maxes(:,s)==0) = 0;
%
%
%
%
%
%         networkconnections(:,s) = thissub_networkconnections;
%    end


%-----------------------------
%Write results

%     out_template.data = zeros(nverts,length(subjects));
%
%     out_template.data(1:size(networkconnections,1),:) = networkconnections;
%     ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh)],out_template);
%
%     out_template.data(1:size(networkconnections,1),:) = dice_diffs;
%     ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicediffs'],out_template);
%
%     out_template.data(1:size(networkconnections,1),:) = dice_maxes;
%     ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicemaxes'],out_template);
%
%     disp(' ')










