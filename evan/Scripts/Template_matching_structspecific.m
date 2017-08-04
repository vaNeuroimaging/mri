% surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
% [subjects, ciftifiles] = textread(surfdatafile,'%s %s');
% 
% tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
% [subjects, tmasks] = textread(tmaskfile,'%s %s');
% 
% subjects = subjects(subnums);
% tmasks = tmasks(subnums);
% ciftifiles = ciftifiles(subnums);

tmaskfile = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/COHORTSELECT_CONTROLS/NEW_TMASKLIST.txt';
[subjects, tmasks] = textread(tmaskfile,'%s %s');
ciftifiles{1} = '/data/cn4/evan/Evan_brain/cifti_timeseries_normalwall/vc39555_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii';

subnums = [1:length(subjects)];%[121:228];

xdistance = 30;

minclustersizemm2 = 30;
minclustersizemm3 = 100;

kden_threshs = .05;


%-----------------------------


%load('/data/cn4/evan/RestingState/Ind_variability/Templates_Yeo.mat')
load('/data/cn4/evan/RestingState/Ind_variability/Templates_consensus.mat')
templates{1}(:,11:12) = [];
IDs{1}(11:12) = [];
templates{1}(:,end-1:end) = [];
IDs{1}(end-1:end) = [];

templates = templates{1};
IDs = IDs{1};

% distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');
% 
% threshdistance = distances > xdistance;
% 
% clear distances

%%
prevstring = [];
for kden_thresh = kden_threshs;
    %disp(kden_thresh)
    
    
    for s = 1:length(subjects)
        
        string = ['Subject ' num2str(subnums(s)) ': loading data'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        %-----------------------------
        %Load data
        
        subject = subjects{s};
        cifti_file = ciftifiles{s};
        
        subdata = ft_read_cifti_mod(cifti_file); 
        
%        [templates, subdata] = register_ciftis_differentvolumes(templates,subdata);
        
        if s==1;
            nverts = size(subdata.data,1);
            out_template = subdata; out_template.data = [];
            out_template.dimord = 'scalar_pos';
            out_template.mapname = cell(1,length(subjects));
            neighbors = cifti_neighbors(out_template);
            dice_diffs = zeros(nverts,length(subjects));
            dice_maxes = zeros(nverts,length(subjects));
            networkconnections = zeros(nverts,length(subjects));
            ncortverts = nnz(out_template.brainstructure==1) + nnz(out_template.brainstructure==2);
            voxvol = out_template.transform(1,1) * out_template.transform(2,2) * out_template.transform(3,3);
            
            
            structinds = ones(nverts,1); structinds(1:ncortverts) = 0;
            structs = unique(structinds);
            
            %threshold templates
            ThreshTemplates = false(size(templates));
            for struct = structs(:)'
                smalltemplates = templates(structinds==struct,:);
                values_sorted = sort(smalltemplates(:),'descend');
                threshval = values_sorted(round(numel(smalltemplates) .* kden_thresh));
                ThreshTemplates(structinds==struct,:) = smalltemplates >= threshval;
            end
            clear values_sorted smalltemplates templates
            
            %get structural divisions
            divisioncounter = 0;
            struct_combination_inds = ones(nverts,'uint8');
            for structnum1 = 1:length(structs)
                for structnum2 = structnum1:length(structs)
                    divisioncounter = divisioncounter +1;
                    struct_combination_inds(structinds==structs(structnum1),structinds==structs(structnum2)) = divisioncounter;
                end
            end
            struct_combination_inds = struct_combination_inds(triu(true(nverts),1));
            tic
        end
        
        out_template.mapname{s} = ['Subject ' num2str(subnums(s))];
        subdata = subdata.data;
        tmask = load(tmasks{s});
        subdata = subdata(:,logical(tmask));
        
        
        %-----------------------------
        %Calculate and threshold correlations
        
        string = ['Subject ' num2str(subnums(s)) ': calculating and thresholding correlation maps'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        
        correlmaps = paircorr_mod(subdata');
        clear subdata
        correlmaps(isnan(correlmaps)) = 0;
        
        matrix_vec=correlmaps(triu(true(nverts),1));
        clear correlmaps
        
        for division = 1:divisioncounter
            
            smallinds = find(struct_combination_inds==division);
            
            %get the matrix values for this pairing of structures
            smallvals = matrix_vec(smallinds);
            
            %sort those values
            [~, smallorder]=sort(smallvals,'descend');
            clear smallvals
            
            %put those values back into the matrix as percentiles within this pairing of structures
            matrix_vec(smallinds(smallorder)) = [true(ceil(numel(smallorder) * kden_thresh),1); false(floor(numel(smallorder) * (1-kden_thresh)),1)];
            clear smallorder smallinds
            
        end
        matrix_vec = logical(matrix_vec);
        
        correlmaps_thresh = false(nverts);
        correlmaps_thresh(triu(true(nverts),1)) = matrix_vec;
        correlmaps_thresh = correlmaps_thresh | (correlmaps_thresh');
        
        clear matrix_vec
        
        
        
        %-----------------------------
        %Match to templates
        
        
        
        dice_coeffs = zeros(nverts,size(ThreshTemplates,2));
        
        for templatenum = 1:size(ThreshTemplates,2);
            string = ['Subject ' num2str(s) ': comparing to network ' num2str(IDs(templatenum))];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            template_big = repmat(ThreshTemplates(:,templatenum)',nverts,1);
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
        
        
        
        string = ['Subject ' num2str(s) ': eliminating patches smaller than ' num2str(minclustersizemm2) 'mm2 or ' num2str(minclustersizemm3) 'mm3'];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        
        hems = {'L','R'};
        for hemnum = 1:length(hems)
            surfaceareafiles{hemnum} = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
            if ~exist(surfaceareafiles{hemnum})
                surffile = [surfaceareafiles{hemnum}(1:end-22) '.surf.gii'];
                system(['wb_command -surface-vertex-areas ' surffile ' ' surfaceareafiles{hemnum}])
            end
            surfaceareas{hemnum} = gifti(surfaceareafiles{hemnum}); 
            surfaceareas{hemnum} = surfaceareas{hemnum}.cdata;
            surfaceareas{hemnum}(out_template.brainstructure((1:length(surfaceareas{hemnum})) + (length(surfaceareas{hemnum}) * (hemnum-1)))==-1) = [];
        end
        
        surfacearea_voxvol = [surfaceareas{1} ; surfaceareas{2} ; (ones(size(correlmaps_thresh,1) - ncortverts , 1) * voxvol)];
        
        temp = zeros(nverts,1);
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
            newblankinds = find(thissub_networkconnections==0);
            if nnz(blankinds)==nnz(newblankinds)
                break
            else
                blankinds = newblankinds;
            end
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









