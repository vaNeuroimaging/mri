surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects, ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects, tmasks] = textread(tmaskfile,'%s %s');

subjects = subjects(1:120);%(121:end);
tmasks = tmasks(1:120);%(121:end);
ciftifiles = ciftifiles(1:120);%(121:end);

outfolder = '/data/cn4/evan/RestingState/Ind_variability/120_infomaps_0002_to_003_in_0002_xd30_BI/';

thresholdarray = [.002: .002: .03];

xdistance = 30;

dmatname = '/data/hcp-zfs/home/laumannt/120_parcellation/modified_cifti_network/normalwall_distmat/distmat_surf_geodesic_vol_euc.mat';

sizethresh = 400;

%% Run Infomap
for s = 1:length(subjects)
    
    disp(['Subject ' num2str(s)])
    subject = subjects{s};
    cifti_file = ciftifiles{s};
    
    subdata = cifti_read(cifti_file);
    tmask = load(tmasks{s});
    subdata = subdata(:,logical(tmask));
    
    rmat = paircorr_mod(subdata');
    clear subdata
    rmat(isnan(rmat)) = 0;
    rmat = FisherTransform(rmat);
    
    suboutfolder = [outfolder subject '/'];
    
    Infomap_new(rmat,dmatname,xdistance,thresholdarray,1,'kden',suboutfolder);
    
    clear rmat
    
    copyfile([suboutfolder '/rawassn.txt'],[outfolder '/' subject '_rawassn.txt'])
    
    rmdir(suboutfolder,'s')
    
end

%% Size-threshold, regularize, and create consensus

%load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
%distances = distances(1:59412,1:59412);
cd(outfolder)

out = zeros(66697,length(subjects));
for s = 1:length(subjects)
    disp(['Subject ' num2str(s)])
    subject = subjects{s};
    simple_assigns = modify_clrfile('simplify',[subject '_rawassn.txt'],sizethresh);
    regularized = rawoutput2clr(simple_assigns);
    regularizedfilename = [subject '_rawassn_minsize' num2str(sizethresh) '_regularize.txt'];
    dlmwrite(regularizedfilename,regularized,'\t');
    out(:,s) = Consensus_networkIDs_fillin2(regularizedfilename,distances);
    %out(:,s) = cifti_read(['Consensus_Powercolors_cleaned_' subject '_rawassn_minsize' num2str(sizethresh) '_regularize.dtseries.nii']);
end

cifti_write_wHDR(out,[],['Consensus_Powercolors_cleaned_allsubs_rawassn_minsize' num2str(sizethresh) '_regularize.dtseries.nii'])