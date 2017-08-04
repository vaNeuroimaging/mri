subjects = {'MSC01','MSC02','MSC03','MSC04'};

%MSCname = 'MSC01';



% system(['wb_command -surface-vertex-areas /data/nil-bluearc/GMT/Laumann/MSC/fs5.3/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/' MSCname '.L.midthickness.32k_fs_LR.surf.gii ' surfaceareafiles{1}])
% system(['wb_command -surface-vertex-areas /data/nil-bluearc/GMT/Laumann/MSC/fs5.3/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/' MSCname '.R.midthickness.32k_fs_LR.surf.gii ' surfaceareafiles{2}])


medial_wall{1} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);
medial_wall{1} = medial_wall{1}.cdata;
medial_wall{2} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);
medial_wall{2} = medial_wall{2}.cdata;

ncortverts = nnz(~medial_wall{1}) + nnz(~medial_wall{2});






%tmaskfile = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/COHORTSELECT_CONTROLS/NEW_TMASKLIST.txt';
%[subjects, tmasks] = textread(tmaskfile,'%s %s');
%ciftifiles{1} = '/data/cn4/evan/Evan_brain/cifti_timeseries_normalwall/vc39555_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii';



xdistance = 30;

minclustersizemm2 = 30;

kden_thresh = .05;

%ncortverts = 59412;

%-----------------------------


%load('/data/cn4/evan/RestingState/Ind_variability/Templates_Yeo.mat')
load('/data/cn4/evan/RestingState/Ind_variability/Templates_consensus.mat')
% templates{1}(:,11:12) = [];
% IDs{1}(11:12) = [];
%templates{1}(:,end-1:end) = [];
%IDs{1}(end-1:end) = [];
% 
% templates = templates{1};
% IDs = IDs{1};

load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat
%%distances = distances(1:ncortverts,1:ncortverts);

threshdistance = distances > xdistance;

clear distances

%%

values_sorted = sort(templates(:),'descend');
threshval = values_sorted(round(numel(templates) .* kden_thresh));
ThreshTemplates = templates >= threshval;
clear values_sorted





prevstring = [];
for s = 1:length(subjects)
    MSCname = subjects{s};
    subject = MSCname;
    %-----------------------------
    %Load data
    
    string = ['Subject ' num2str(s) ': calculating correlation maps'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    
    tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
    [sessions tmasks] = textread(tmaskfile,'%s %s');
    for i = 1:length(sessions)
        if i == 1
            subdata = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' sessions{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
            
            tmask = load(tmasks{i});
            
            if s==1;
                out_template = subdata; out_template.data = [];
                out_template.dimord = 'scalar_pos';
                out_template.mapname = cell(1,length(subjects));
                
                nverts = size(subdata.data,1);
                ncortverts = nnz(out_template.brainstructure==1) + nnz(out_template.brainstructure==2);
                neighbors = cifti_neighbors(out_template);
                neighbors = neighbors(1:ncortverts,:);
                threshdistance = threshdistance(1:ncortverts,1:ncortverts);
                ThreshTemplates = ThreshTemplates(1:ncortverts,:);
                dice_diffs = zeros(ncortverts,length(subjects));
                dice_maxes = zeros(ncortverts,length(subjects));
                networkconnections = zeros(ncortverts,length(subjects));
                voxvol = out_template.transform(1,1) * out_template.transform(2,2) * out_template.transform(3,3);
            end
            subdata = subdata.data;
            
        else
            temp = ft_read_cifti_mod(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' sessions{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
            subdata = [subdata temp.data];
            clear temp
            tmask = [tmask; load(tmasks{i})];
        end
    end
    subdata = subdata(1:ncortverts,logical(tmask));
    
    
    out_template.mapname{s} = MSCname;
    
    
    
    %-----------------------------
    %Calculate and threshold correlations
    
    
    correlmaps = paircorr_mod(subdata');
    clear subdata
    correlmaps(isnan(correlmaps)) = 0;
    
    temp=correlmaps(triu(true(ncortverts),1));
    v=sort(temp,'descend');
    clear temp
    r_sub_thresh = v(round(kden_thresh * numel(v)));
    clear v
    correlmaps_thresh = correlmaps > r_sub_thresh;
    
    clear correlmaps
    
    
    %-----------------------------
    %Match to templates
    
    
    
    dice_coeffs = zeros(ncortverts,size(ThreshTemplates,2));
    
    for templatenum = 1:size(ThreshTemplates,2);
        string = ['Subject ' num2str(s) ': comparing to network ' num2str(IDs(templatenum))];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        template_big = repmat(ThreshTemplates(:,templatenum)',ncortverts,1);
        dice_coeffs(:,templatenum) = sum((template_big & correlmaps_thresh & threshdistance),2) ./ sum(((template_big | correlmaps_thresh) & threshdistance),2);
    end
    clear template_big
    
    [sorteddice, sorti] = sort(dice_coeffs,2,'descend');
    sorti = sorti(:,1)';
    dice_maxes(:,s) = sorteddice(:,1);
    dice_diffs(:,s) = sorteddice(:,1) - sorteddice(:,2);
    clear sorteddice
    
    thissub_networkconnections = IDs(sorti)';
    
    thissub_networkconnections(dice_maxes(:,s)==0) = 0;
    
    networkconnections(:,s) = thissub_networkconnections;
    
    clear dice_coeffs
    
    
end

%-----------------------------
%Write results

out_template.data = zeros(nverts,length(subjects));

out_template.data(1:size(networkconnections,1),:) = networkconnections;
ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh)],out_template);

%out_template.data(1:size(networkconnections,1),:) = dice_diffs;
%ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicediffs'],out_template);

%out_template.data(1:size(networkconnections,1),:) = dice_maxes;
%ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicemaxes'],out_template);

disp(' ')
prevstring = [];


%%



%-----------------------------
%Size threshold

for s = 1:length(subjects)
    
    thissub_networkconnections = networkconnections(:,s);
    
    string = ['Subject ' num2str(s) ': eliminating patches smaller than ' num2str(minclustersizemm2) 'mm2'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    
    hems = {'L','R'};
    for hemnum = 1:length(hems)
        %surfaceareafiles{hemnum} = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
        %surfaceareafiles{hemnum} = ['/data/hcp-zfs/shared-nil/adeyemob/HCP_MAIN/OTHER/fsaverage_LR32K/' subject '/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
        %surfaceareafiles{hemnum} = ['/data/nil-bluearc/GMT/Evan/CoE/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
        surfaceareafiles{hemnum} = ['/data/nil-bluearc/GMT/Laumann/MSC/fs5.3_native/FREESURFER_fs_LR/' MSCname '/7112b_fs_LR/fsaverage_LR32k/' MSCname '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
        
        if ~exist(surfaceareafiles{hemnum})
            surffile = [surfaceareafiles{hemnum}(1:end-22) '.surf.gii'];
            system(['wb_command -surface-vertex-areas ' surffile ' ' surfaceareafiles{hemnum}])
        end
        surfaceareas{hemnum} = gifti(surfaceareafiles{hemnum});
        surfaceareas{hemnum} = surfaceareas{hemnum}.cdata;
        surfaceareas{hemnum}(out_template.brainstructure((1:length(surfaceareas{hemnum})) + (length(surfaceareas{hemnum}) * (hemnum-1)))==-1) = [];
    end
    
    surfacearea_voxvol = [surfaceareas{1} ; surfaceareas{2} ; (ones(size(correlmaps_thresh,1) - ncortverts , 1) * voxvol)];
    
    temp = zeros(ncortverts,1);
    for ID = IDs(:)'
        clustereddata = cifti_cluster_surfacearea_volume(thissub_networkconnections,ID-.5,ID+.5,minclustersizemm2,minclustersizemm2,ncortverts,surfacearea_voxvol,neighbors);
        clustereddata(logical(clustereddata)) = ID;
        temp = temp + clustereddata;
    end
    
    thissub_networkconnections = temp;
    
    
    
    %-----------------------------
    %Fill in removed spots
    
    
    
    blankinds = find(thissub_networkconnections==0);
    while ~isempty(blankinds)
        temp = thissub_networkconnections;
        for ind = blankinds(:)'
            indneighs = neighbors(ind,2:end); indneighs(isnan(indneighs)) = [];
            neighvals = thissub_networkconnections(indneighs); neighvals(neighvals==0) = [];
            if ~isempty(neighvals)
                temp(ind) = mode(neighvals);
            end
        end
        thissub_networkconnections = temp;
        blankinds = find(thissub_networkconnections==0);
    end
    
    thissub_networkconnections(dice_maxes(:,s)==0) = 0;
    
    
    
    networkconnections(:,s) = thissub_networkconnections;
end

out_template.data(1:size(networkconnections,1),:) = networkconnections;
ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_minsize' num2str(minclustersizemm2) 'mm2'],out_template);















