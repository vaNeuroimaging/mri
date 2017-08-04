function Template_matching_cortexonly_parcels_onesub(subdata,parcels,distances,xdistance,kden_thresh,templatesfile,outname)
%Template_matching_cortexonly_onesub(subdata,parcels,distances,xdistance,kden_thresh,templatesfile,outname)


if ischar(subdata)
    subdata = ft_read_cifti_mod(subdata);
end

if ischar(distances)
    distances = smartload(distances);
    distances = uint8(distances);
end

if ischar(parcels)
    parcels = ft_read_cifti_mod(parcels); parcels = parcels.data;
end


%% Parcels
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];
parcel_centroids = zeros(length(parcelIDs),1);
parcel_tcs = zeros(size(subdata.data,2),length(parcelIDs));
for p = 1:length(parcelIDs)
    parcelinds = find(parcels==parcelIDs(p));
    [~,centroidind] = min(sum(distances(parcelinds,parcelinds),2));
    parcel_centroids(p) = parcelinds(centroidind);
    
    parcel_tcs(:,p) = mean(subdata.data(parcelinds,:),1);
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
threshdistance = threshdistance(1:ncortverts,parcel_centroids);
ThreshTemplates = ThreshTemplates(1:ncortverts,:);


subdata = subdata.data(1:ncortverts,:);


%-----------------------------
%Calculate and threshold vertexwise correlations

correlmaps = paircorr_mod(subdata',parcel_tcs);
correlmaps(isnan(correlmaps)) = 0;
v=sort(correlmaps(:),'descend');
r_sub_thresh = v(round(kden_thresh * numel(v)));
clear v
correlmaps_thresh = correlmaps > r_sub_thresh;




% %correlmaps = paircorr_mod(subdata');
% clear subdata
% correlmaps(isnan(correlmaps)) = 0;
% 
% temp=correlmaps(triu(true(ncortverts),1));
% v=sort(temp,'descend');
% clear temp
% r_sub_thresh = v(round(kden_thresh * numel(v)));
% clear v
% correlmaps_thresh = correlmaps > r_sub_thresh;

clear correlmaps


%-----------------------------
%Match to templates



dice_coeffs = zeros(length(parcelIDs),size(ThreshTemplates,2));

for templatenum = 1:size(ThreshTemplates,2);
    string = ['Comparing to network ' num2str(IDs(templatenum))];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    template_big = repmat(ThreshTemplates(:,templatenum),1,length(parcelIDs));
    dice_coeffs(:,templatenum) = sum((template_big & correlmaps_thresh & threshdistance),1) ./ sum(((template_big | correlmaps_thresh) & threshdistance),1);
end
clear template_big

[sorteddice, sorti] = sort(dice_coeffs,2,'descend');
sorti = sorti(:,1)';
dice_maxes = sorteddice(:,1);
dice_diffs = sorteddice(:,1) - sorteddice(:,2);
%clear sorteddice

thissub_networkconnections = IDs(sorti)';

thissub_networkconnections(dice_maxes==0) = 0;

%networkconnections = thissub_networkconnections;

%clear dice_coeffs



out_template.data = zeros(size(parcels,1),1);
out_template_dicemaxes = out_template;
out_template_dicevals = out_template; out_template_dicevals.data = zeros(size(parcels,1),size(dice_coeffs,2));
for p = 1:length(parcelIDs)
    out_template.data(parcels==parcelIDs(p)) = thissub_networkconnections(p);
    out_template_dicemaxes.data(parcels==parcelIDs(p)) = dice_maxes(p);
    out_template_dicevals.data(parcels==parcelIDs(p),:) = repmat(dice_coeffs(p,:),nnz(parcels==parcelIDs(p)),1);
end


%-----------------------------
%Write results

%out_template.data = zeros(nverts,1);

%out_template.data(1:size(networkconnections,1),:) = networkconnections;
ft_write_cifti_mod([outname '_templatematch_dice_kden' num2str(kden_thresh) '_parcels'],out_template);

%out_template.data(1:size(networkconnections,1),:) = dice_diffs;
%ft_write_cifti_mod(['Templatematch_dice_bysubject_kden' num2str(kden_thresh) '_dicediffs'],out_template);

%out_template.data(1:size(networkconnections,1),:) = dice_maxes;
ft_write_cifti_mod([outname '_templatematch_dice_bysubject_kden' num2str(kden_thresh) '_parcels_dicemaxes'],out_template_dicemaxes);

ft_write_cifti_mod([outname '_templatematch_dice_bysubject_kden' num2str(kden_thresh) '_parcels_dicevals'],out_template_dicevals);


disp(' ')


















