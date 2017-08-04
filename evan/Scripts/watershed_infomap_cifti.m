parcelfile = '/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii';

corrmatname = 'corrmat.mat';
outputfolder = ['/data/cn4/evan/Temp/partials/partial_corr/'];

thresholdarray = [.005 : .001 : .05];
xdistance = 30;
networksizeminimum = 5;

outputstem = 'Partials';

calc_corrmat = 0;
calc_distances = 1;


cohortfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects surfdatafile] = textread(cohortfile,'%s %s');
tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

distancesfile = '/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat';


 
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


%% Parcel networks

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


%% Find watershed parcel centers

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
    
    save([outputfolder '/parcel_distances.mat'],'parcel_distances')
    
else
    parcel_distances = smartload([outputfolder '/parcel_distances.mat']);
end


%% Run Infomap

Infomap_new(all_parcel_corrmat, [outputfolder '/parcel_distances.mat'], xdistance, thresholdarray, 0, 'kden', outputfolder)


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





regularized = rawoutput2clr(simple_assigns);
regularized(regularized < 2) = 0;
regularized = regularized -1;
dlmwrite(['rawassn_minize' num2str(networksizeminimum) '_regularized.txt'],regularized);

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





    
    
 