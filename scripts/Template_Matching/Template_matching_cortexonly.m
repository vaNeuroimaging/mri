function Template_matching_cortexonly(template_matching_params_file)
%Template_matching_cortexonly(template_matching_params_file)
%
% Creates subject versions of template systems by assessing the similarity
% of each vertex's connectivity map (as the dice overlap of thresholded
% binarized maps) in each subject to a set of template system connectivity
% maps. Each vertex is assigned the identity of the system with the most
% similar template connectivity map. Subject maps are then size-thresholded
% to eliminate contiguous patches smaller than a given size (in mm2); these
% removed spots are then filled in with neighboring systems.
%
% Requires a parameters file (a .m file) which will be executed to load
% needed parameters, including:
%
% a datalist and tmasklist for the subjects to be run
% a geodesic exclusion distance (connectivity map similarity is not
%  assessed within this distance of the seed vertex)
% a minimum size below which contiguous system patches will be eliminated
% a density threshold to binarize at
% a set of template system connectivity maps
% a point-to-point distances matrix
%
%
% Requires the Cifti Resources scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ and subfolders)
%
%EMG 06/26/15


%Find parameters file
[paramspath,paramsname,paramsextension] = fileparts(template_matching_params_file);
origpath = pwd;
if ~isempty(paramspath)
    cd(paramspath)
end

%Load parameters
params = feval(paramsname);
varnames = fieldnames(params);
for i = 1:length(varnames)
    evalc([varnames{i} ' = params.' varnames{i}]);
end
clear varnames params


%% Load templates, distances, and subject lists

[subjects, ciftifiles] = textread(surfdatafile,'%s %s');

load(templatesfile)

distances = smartload(distancesfile);

threshdistance = distances > xdistance;

clear distances

values_sorted = sort(templates(:),'descend');
threshval = values_sorted(round(numel(templates) .* kden_thresh));
ThreshTemplates = templates >= threshval;
clear values_sorted


%% Loop through subjects
prevstring = [];
for s = 1:length(subjects)
    
    %-----------------------------
    %Load data
    
    string = ['Subject ' num2str(s) ': thresholding correlation maps'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    subject = subjects{s};
    cifti_file = ciftifiles{s};
    
    correlmaps = ft_read_cifti_mod(cifti_file);
    
    
    %-----------------------------
    %Set up needed variables
    if s==1;
        out_template = correlmaps; out_template.data = [];
        out_template.dimord = 'scalar_pos';
        out_template.mapname = cell(1,length(subjects));
        
        nverts = size(correlmaps.data,1);
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
    
    out_template.mapname{s} = subject;
    correlmaps = correlmaps.data;
    correlmaps = correlmaps(1:ncortverts,1:ncortverts);
    
    
    %-----------------------------
    %Calculate and threshold vertexwise correlations
    
    
    %correlmaps = paircorr_mod(correlmaps');
    %clear subdata
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

out_template.data(1:size(networkconnections,1),:) = dice_diffs;
ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicediffs'],out_template);

out_template.data(1:size(networkconnections,1),:) = dice_maxes;
ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicemaxes'],out_template);

disp(' ')
prevstring = [];






%-----------------------------
%Size threshold results

for s = 1:length(subjects)
    
    thissub_networkconnections = networkconnections(:,s);
    
    string = ['Subject ' num2str(s) ': eliminating patches smaller than ' num2str(minclustersizemm2) 'mm2'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    %Find / make vertexwise surface area files
    hems = {'L','R'};
    for hemnum = 1:length(hems)
        surfaceareafiles{hemnum} = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.' hems{hemnum} '.midthickness.32k_fs_LR_surfaceareas.func.gii'];
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


disp(' ')












