function Template_matching_cortexonly_onesub(subdata,distances,xdistance,kden_thresh,templatesfile,outname)
%Template_matching_cortexonly_onesub(subdata,distances,xdistance,kden_thresh,templatesfile,outname)


if ischar(subdata)
    subdata = ft_read_cifti_mod(subdata);
end

if ischar(distances)
    distances = smartload(distances);
    distances = uint8(distances);
end


%% Load templates, distances, and subject lists



load(templatesfile)


threshdistance = distances > xdistance;

clear distances

values_sorted = sort(templates(:),'descend');
threshval = values_sorted(round(numel(templates) .* kden_thresh));
ThreshTemplates = templates >= threshval;
clear values_sorted


%% Loop through subjects
prevstring = [];


%-----------------------------
%Load data

string = ['Calculating correlation maps'];
fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
prevstring = string;



out_template = subdata; out_template.data = [];

nverts = size(subdata.data,1);
ncortverts = nnz(out_template.brainstructure==1) + nnz(out_template.brainstructure==2);
threshdistance = threshdistance(1:ncortverts,1:ncortverts);
ThreshTemplates = ThreshTemplates(1:ncortverts,:);


subdata = subdata.data(1:ncortverts,:);


%-----------------------------
%Calculate and threshold vertexwise correlations


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
    string = ['Comparing to network ' num2str(IDs(templatenum))];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    template_big = repmat(ThreshTemplates(:,templatenum)',ncortverts,1);
    dice_coeffs(:,templatenum) = sum((template_big & correlmaps_thresh & threshdistance),2) ./ sum(((template_big | correlmaps_thresh) & threshdistance),2);
end
clear template_big

[sorteddice, sorti] = sort(dice_coeffs,2,'descend');
sorti = sorti(:,1)';
dice_maxes = sorteddice(:,1);
dice_diffs = sorteddice(:,1) - sorteddice(:,2);
clear sorteddice

thissub_networkconnections = IDs(sorti)';

thissub_networkconnections(dice_maxes==0) = 0;

networkconnections = thissub_networkconnections;

clear dice_coeffs



%-----------------------------
%Write results

out_template.data = zeros(nverts,1);

out_template.data(1:size(networkconnections,1),:) = networkconnections;
ft_write_cifti_mod([outname '_templatematch_jaccard_kden' num2str(kden_thresh)],out_template);
set_cifti_powercolors([outname '_templatematch_jaccard_kden' num2str(kden_thresh) '.dtseries.nii'])

out_template.data(1:size(networkconnections,1),:) = dice_diffs;
ft_write_cifti_mod([outname '_templatematch_dice_bysubject_kden' num2str(kden_thresh) '_jaccarddiffs'],out_template);

out_template.data(1:size(networkconnections,1),:) = dice_maxes;
ft_write_cifti_mod([outname '_templatematch_dice_bysubject_kden' num2str(kden_thresh) '_jaccardmaxes'],out_template);


disp(' ')


















