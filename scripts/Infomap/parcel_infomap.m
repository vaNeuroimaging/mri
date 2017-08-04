function parcel_infomap(parcel_infomap_params_file)
%parcel_infomap(parcel_infomap_params_file)
%
% Runs infomap on a set of parcels at a variety of thresholds,
% size-thresholds and regularizes the resulting communities, and writes the
% communities back into cifti format for easy viewing.
%
% Requires a parameters file (a .m file) which will be executed to load
% needed parameters, including:
%
% the parcels being evaluated
% the output folder
% thresholds to run
% geodesic distance exclusion
% minimum community size
% whether or not to calculate new parcel-to-parcel correlations from cifti
%  timecourses, and the files needed to do that 
% whether or not to calculate new parcel-to-parcel geodesic distances, and
%  the files needed to do that 
%
%
% Requires the Cifti Resources scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ and subfolders)
%
%
%EMG 06/25/15

%% Load parameters
%Find params file
[paramspath,paramsname,paramsextension] = fileparts(parcel_infomap_params_file);
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


cd(origpath)


%% Load files

 
[subjects surfdatafile] = textread(cohortfile,'%s %s');
[subjects tmasks] = textread(tmaskfile,'%s %s');


mkdir(outputfolder)

parcels = ft_read_cifti_mod(parcelfile);
cifti_template = parcels; cifti_template.data = [];
cifti_template.dimord = 'scalar_pos';
cifti_template.mapname = cell(1,length(thresholdarray));
for col = 1:length(thresholdarray)
    cifti_template.mapname{col} = ['kden ' num2str(thresholdarray(col)) ', xdistance ' num2str(xdistance)];
end


parcels = parcels.data;
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];


%% Calculate the parcel correlation matrix

if calc_corrmat
    
    parcel_corrmat = zeros(length(parcelIDs),length(parcelIDs),length(subjects));
    
    for s = 1:length(subjects)
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        
        tmask = load(tmasks{s});
        
        
        surf_timecourses = ft_read_cifti_mod(surfdatafile{s});
        
        
        surf_timecourses = surf_timecourses.data(:,logical(tmask));
        surf_timecourses(isnan(surf_timecourses)) = 0;
        
        parceltimecourses = zeros(length(parcelIDs),size(surf_timecourses,2));
        for IDnum = 1:length(parcelIDs)
            ID = parcelIDs(IDnum);
            parceltimecourses(IDnum,:) = nanmean(surf_timecourses(parcels==ID,:),1);
        end
        
        parcel_corrmat(:,:,s) = FisherTransform(paircorr_mod(parceltimecourses'));
                
    end
    
    all_parcel_corrmat = nanmean(parcel_corrmat,3);
    
    save([outputfolder '/' corrmatname],'all_parcel_corrmat')
    
else
    all_parcel_corrmat = smartload([outputfolder '/' corrmatname]);
end


%% Find parcel centers

if calc_distances
    
    
    distances = smartload(distancesfile);
    distances = uint8(distances(1:length(parcels),1:length(parcels)));
    
    centers = zeros(length(parcelIDs),1);
    
    for IDnum = 1:length(parcelIDs)
        ID = parcelIDs(IDnum);
        parcelverts = find(parcels==ID);
        thisparcel_distances = sum(distances(parcelverts,parcelverts),2);
        [ign,mini] = min(thisparcel_distances);
        centers(IDnum) = parcelverts(mini);
    end
    
    parcel_distances = distances(centers,centers);
    
    
    clear distances
    
    save([outputfolder '/' parceldistancesname],'parcel_distances')
    
else
    parcel_distances = smartload([outputfolder '/' parceldistancesname]);
end


%% Run Infomap

Run_Infomap(all_parcel_corrmat, [outputfolder '/' parceldistancesname], xdistance, thresholdarray, 0, outputfolder)


%% Modify color assignments with miminum network size criterion, regularize, and write out results


cd(outputfolder)

rawassns = load('rawassn.txt');
rawassns_onbrain = zeros(size(parcels,1),size(rawassns,2));
for IDnum = 1:length(parcelIDs)
        ID = parcelIDs(IDnum);
        parcelverts = find(parcels==ID);
        for col = 1:size(rawassns,2)
            rawassns_onbrain(parcelverts,col) = rawassns(IDnum,col);
        end
end
cifti_template.data = rawassns_onbrain;
ft_write_cifti_mod(['rawassn'],cifti_template);



simple_assigns = modify_clrfile('simplify','rawassn.txt',networksizeminimum);

simple_assigns_onbrain = zeros(size(parcels,1),size(simple_assigns,2));
for IDnum = 1:length(parcelIDs)
        ID = parcelIDs(IDnum);
        parcelverts = find(parcels==ID);
        for col = 1:size(simple_assigns,2)
            simple_assigns_onbrain(parcelverts,col) = simple_assigns(IDnum,col);
        end
end
cifti_template.data = simple_assigns_onbrain;
ft_write_cifti_mod(['rawassn_minsize' num2str(networksizeminimum)],cifti_template);





regularized = regularize(simple_assigns);
dlmwrite(['rawassn_minsize' num2str(networksizeminimum) '_regularized.txt'],regularized);

regularized_onbrain = zeros(size(parcels,1),size(simple_assigns,2));
for IDnum = 1:length(parcelIDs)
        ID = parcelIDs(IDnum);
        parcelverts = find(parcels==ID);
        for col = 1:size(simple_assigns,2)
            regularized_onbrain(parcelverts,col) = regularized(IDnum,col);
        end
end
cifti_template.data = regularized_onbrain;
ft_write_cifti_mod(['rawassn_minsize' num2str(networksizeminimum) '_regularized'],cifti_template);





    
    
 