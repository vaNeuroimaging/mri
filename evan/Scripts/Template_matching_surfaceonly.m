surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

subnums = [1:length(subjects)];%[121:228];

subjects = subjects(subnums);
tmasks = tmasks(subnums);
ciftifiles = ciftifiles(subnums);

xdistance = 30;

minclustersizemm2 = 30;

kden_threshs = .05;


%-----------------------------


load('/data/cn4/evan/RestingState/Ind_variability/Templates_Yeo.mat')
%load('/data/cn4/evan/RestingState/Ind_variability/Templates_consensus.mat')
% templates{1}(:,11:12) = [];
% IDs{1}(11:12) = [];
% templates{1}(:,end-1:end) = [];
% IDs{1}(end-1:end) = [];


load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
 
threshdistance = distances > xdistance;
 
clear distances

%%

for kden_thresh = kden_threshs;
    %disp(kden_thresh)
    
    
    
    
    
    prevstring = [];
    for s = 1:length(subjects)
        
        
        %-----------------------------
        %Load data
        
        string = ['Subject ' num2str(subnums(s)) ': calculating correlation maps'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        subject = subjects{s};
        cifti_file = ciftifiles{s};
        
        subdata = ft_read_cifti_mod(cifti_file); 
        if s==1;
            nverts = size(subdata.data,1);
            out_template = subdata; out_template.data = [];
            out_template.dimord = 'scalar_pos';
            out_template.mapname = cell(1,length(subjects));
            ncortverts = nnz(out_template.brainstructure==1) + nnz(out_template.brainstructure==2);
            networkconnections = zeros(ncortverts,length(subjects));
            neighbors = cifti_neighbors(out_template); neighbors = neighbors(1:ncortverts,:);
            dice_diffs = zeros(ncortverts,length(subjects));
            dice_maxes = zeros(ncortverts,length(subjects));
            
            threshdistance = threshdistance(1:ncortverts,1:ncortverts);
            
            templates = templates{1}(1:ncortverts,:);
            IDs = IDs{1};
            values_sorted = sort(templates(:),'descend');
            threshval = values_sorted(round(numel(templates) .* kden_thresh));
            ThreshTemplates = templates >= threshval;
            clear values_sorted    
            
            divisions = 4;
            for i = 1:(divisions-1)
                div_verts{i} = [(round(ncortverts / divisions) * (i-1)) + 1 : (round(ncortverts / divisions) * i)];
            end
            div_verts{divisions} = [(round(ncortverts / divisions) * (divisions-1)) + 1 : ncortverts];

        end
        out_template.mapname{s} = string;
        subdata = subdata.data;
        tmask = load(tmasks{s});
        subdata = subdata(1:ncortverts,logical(tmask));
        
        
        %-----------------------------
        %Calculate and threshold correlations
        
        
        correlmaps = zeros(ncortverts);
        for i = 1:divisions
            for j = i:divisions
                mat = paircorr_mod(subdata(div_verts{i},:)',subdata(div_verts{j},:)'); 
                correlmaps(div_verts{i},div_verts{j}) = mat;
                if i~=j
                    correlmaps(div_verts{j},div_verts{i}) = mat';
                end
            end
        end
        clear mat
        
        %correlmaps = paircorr_mod(subdata');
        clear subdata
        correlmaps(isnan(correlmaps)) = 0;
        
        correlmaps=correlmaps(triu(true(ncortverts),1));
        v=sort(correlmaps,'descend');
        r_sub_thresh = v(round(kden_thresh * numel(v)));
        clear v
        correlmaps_thresh_vec = correlmaps > r_sub_thresh;
        clear correlmaps
        correlmaps_thresh = false(ncortverts);
        correlmaps_thresh(triu(true(ncortverts),1)) = correlmaps_thresh_vec;
        clear correlmaps_thresh_vec
        correlmaps_thresh = correlmaps_thresh | correlmaps_thresh';
        
        
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
        
        clear dice_coeffs
        
        
        
        %-----------------------------
        %Size threshold
        
        
        
        string = ['Subject ' num2str(s) ': eliminating patches smaller than ' num2str(minclustersizemm2) 'mm2'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        surfaceareafile = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.LR.midthickness.32k_fs_LR_surfaceareas.dtseries.nii'];
        if ~exist(surfaceareafile)
                    
        hems = {'L','R'};
        for hemnum = 1:length(hems)
            surfaceareafiles{hemnum} = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
            surffile = [surfaceareafiles{hemnum}(1:end-22) '.surf.gii'];
            system(['wb_command -surface-vertex-areas ' surffile ' ' surfaceareafiles{hemnum}])
            system(['wb_command -cifti-create-dense-timeseries /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.LR.midthickness.32k_fs_LR_surfaceareas.dtseries.nii -left-metric ' surfaceareafiles{1} ' -roi-left /data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii -right-metric ' surfaceareafiles{2} ' -roi-right /data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii'])
            
        end
        end
        
        
        surfacearea = ft_read_cifti_mod(surfaceareafile); surfacearea = surfacearea.data;
        
        temp = zeros(ncortverts,1);
        for ID = IDs(:)'
            clustereddata = cifti_cluster_surfacearea_volume(thissub_networkconnections,ID-.5,ID+.5,minclustersizemm2,0,ncortverts,surfacearea,neighbors);
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
            
        
        
        
        
        networkconnections(:,s) = thissub_networkconnections;
    end
    
    
    %-----------------------------
    %Write results
    
    out_template.data = zeros(nverts,length(subjects));
    
    out_template.data(1:size(networkconnections,1),:) = networkconnections;
    ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh)],out_template);
    
    out_template.data(1:size(networkconnections,1),:) = dice_diffs;
    ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicediffs'],out_template);
    
    out_template.data(1:size(networkconnections,1),:) = dice_maxes;
    ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicemaxes'],out_template);
    
    disp(' ')
end









