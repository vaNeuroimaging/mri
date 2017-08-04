%% Subjects

%subjects = '/home/data/subjects/processing_list_10min.txt';
%subjects = '/home/data/subjects/processing_list_121415.txt';
%subjects = '/home/data/subjects/list_forJDP.txt';
%subjects = {'MAV008'};%,'MAV008'};%'MAV014','MAV016','MAV024','MAV033'};
%subjects = '/home/data/subjects/DART_010416_3.txt';
%subjects = '/home/data/subjects/processing_list_041216.txt';
subjects = '/home/data/subjects/DART.txt';
%subjects = {'Movie1'};

%% Steps to run

pairwise_correlations = 1;
template_matching = 0;
infomap = 0;
parcellation = 0;
homogeneity_testing = 0;
parcel_infomap_and_hub_detection = 0;
groupparcel_infomap_and_hub_detection = 0;




%% Parameters

%functional_sequences = {'restingstate'};
functional_sequences = {'RSFC'};%{'Movie'};%

geodesic_smooth = 2.55;%4.25;%

xdist = 30;

infomapthresholds = [.02 : .005 : .05];
infomap_minsize = 400;
            
parcel_infomapthresholds = [.01:.005:.2];
parcel_infomap_minsize = 5;
parcel_infomap_xdist = 30;

parcellation_subsample = 100;

parcel_borderthresh = .5;

homogeneity_testing_iterations = 100;
homogeneitythresh = 60;
groupparcelfile = '/home/data/atlases/Group_parcellation/Parcels_LR.dlabel.nii';
groupparcelcommunitiesfile = '/home/data/atlases/Group_parcellation/Parcel_Communities.dtseries.nii';


templatesfile = '/home/data/scripts/Template_Matching/Templates_consensus.mat';
templates_kden = .05;

hub_detection_kden = .1;


%% Get subjects

if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end



subcount = length(subjects);







disp('ANALYSIS')
    
for subnum = 1:subcount
    
    systemresult = cell(0,2);
    
    subject = subjects{subnum};
    
    warning off
    
    %try
        
        %% Set up folders
        
        ciftifolder = ['/home/data/subjects/' subject '/cifti/']; 
        correlationfolder = ['/home/data/subjects/' subject '/connectome/']; 
        infomapfolder = ['/home/data/subjects/' subject '/infomap/'];
        parcellationfolder = ['/home/data/subjects/' subject '/parcellation/'];
        parcel_infomapfolder = ['/home/data/subjects/' subject '/parcellation/parcel_infomap/'];
        groupparcel_infomapfolder = ['/home/data/subjects/' subject '/groupparcel_infomap/'];
        homogeneityfolder = ['/home/data/subjects/' subject '/parcellation/homogeneity_testing/'];
        templatematchingfolder = ['/home/data/subjects/' subject '/template_matching/'];
        logfolder = ['/home/data/subjects/' subject '/processing_logs/']; 
        
        logfile2 = [logfolder datestr(now) '.txt']; tokens = tokenize(logfile2,' '); logfile2 = [tokens{1} '_' tokens{2}];
        
        
        
        if pairwise_correlations
            parcels = ft_read_cifti_mod(groupparcelfile);
            parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
            communities = ft_read_cifti_mod(groupparcelcommunitiesfile);
            parcels.data = parcels.data(communities.brainstructure(1:min(length(communities.brainstructure),length(parcels.brainstructure))) > 0);
            parcel_communities = zeros(length(parcelIDs),1);
            for IDnum = 1:length(parcelIDs)
                parcel_communities(IDnum,1) = mode(communities.data(parcels.data==parcelIDs(IDnum)));
            end
            
            [~,temp,~] = fileparts(groupparcelfile);
            [~,parcellation_shortname,~] = fileparts(temp);
            
            for seq = 1:length(functional_sequences)
                disp(['Subject ' subject ', ' functional_sequences{seq} ': running group parcel correlations']);
                ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                mkdir([ciftifolder '/parcel_timecourses/'])
                parcel_timecourse_filename = [ciftifolder '/parcel_timecourses/' functional_sequences{seq} '_' parcellation_shortname '.ptseries.nii'];
                system(['wb_command -cifti-parcellate ' ciftifile ' ' groupparcelfile ' COLUMN ' parcel_timecourse_filename]);
                data = ft_read_cifti_mod(parcel_timecourse_filename);
                timecourses = data.data;
                data.dimord = 'chan_chan';
                data.data = paircorr_mod(data.data');
                data.data = FisherTransform(data.data);
                data.data(isnan(data.data)) = 0;
                ft_write_cifti_mod([correlationfolder functional_sequences{seq} '_' parcellation_shortname '_corr'],data);

                
                title_string = [subject ': ' parcellation_shortname ' corrmat'];
                title_string(strfind(title_string,'_')) = ' ';
                parcel_correlmat_figmaker(data.data,parcel_communities,[-.75 .75],title_string);
                
                export_fig(gca,[correlationfolder parcellation_shortname '_corr.pdf'])

            end
            
        end
        
        
        
        %% Template Matching
        
        if template_matching
            mkdir(templatematchingfolder)
            cd(templatematchingfolder)
            dlmwrite(logfile2,'running template matching','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': running template matching']);
                ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                dmatfile = [ciftifolder '/distances/normalwall_distmat_uint8.mat'];
                if ~exist(dmatfile)
                    dmatfile = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_standardsubcort_uint8.mat';
                end
                distances = smartload(dmatfile);
                Template_matching_cortexonly_onesub(ciftifile,distances,xdist,templates_kden,templatesfile,functional_sequences{seq})
                
                %parcelfile = [parcellationfolder '/' functional_sequences{seq} '_parcels_edgethresh_' num2str(parcel_borderthresh) '.dtseries.nii'];
                %Template_matching_cortexonly_parcels_onesub(ciftifile,parcelfile,distances,xdist,templates_kden,templatesfile,functional_sequences{seq})
                
            end
        end
        
        
        
        %% Infomap
        
        if infomap
            
            mkdir(infomapfolder)
            
            dlmwrite(logfile2,'running infomap','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                
                outfolder = [infomapfolder '/' functional_sequences{seq} '/'];
                mkdir(outfolder)
            
                disp(['Subject ' subject ', ' functional_sequences{seq} ': running infomap']);
                
                corrfile = [correlationfolder '/' functional_sequences{seq} '_corr.dconn.nii'];
                
                dmatfile = [ciftifolder '/distances/normalwall_distmat.dconn.nii'];
                if ~exist(dmatfile)
                    dmatfile = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_standardsubcort_uint8.mat';
                end
                
                Run_Infomap(corrfile, dmatfile, xdist, infomapthresholds, 0, outfolder);
                
                
                simplified=modify_clrfile('simplify',[outfolder '/rawassn.txt'],infomap_minsize);
                regularized =  regularize(simplified,[outfolder '/rawassn_minsize' num2str(infomap_minsize) '_regularized.txt']);
                
                outcifti = ft_read_cifti_mod([ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii']);
                outcifti.data = regularized;
                outcifti.dimord = 'scalar_pos';
                for i = 1:length(infomapthresholds)
                    outcifti.mapname{1,i} = ['minsize=' num2str(infomap_minsize) '; xdist=' num2str(xdist) ';thresh=' num2str(infomapthresholds(i))];
                end
                ft_write_cifti_mod([outfolder '/rawassn_minsize' num2str(infomap_minsize) '_regularized.dscalar.nii'],outcifti);
                
                consensus_maker_knowncolors([outfolder '/rawassn_minsize' num2str(infomap_minsize) '_regularized.dscalar.nii'])
                
                
            end
        end
        
        %% Parcellation
        
        if parcellation
            
            mkdir(parcellationfolder)
            cd(parcellationfolder)
            
            dlmwrite(logfile2,'running parcellation','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': running parcellation']);
                
                mkdir(['/home/data/subjects/' subject '/connectome/'])
                ciftifile = [ciftifolder '/cifti_timeseries_smallwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                data = ft_read_cifti_mod(ciftifile);
                
                corr = paircorr_mod(data.data');
                corr = FisherTransform(corr);
                corr(isnan(corr)) = 0;
                
                data.dimord = 'pos_pos';
                data.data = corr;
                clear corr
                
                ft_write_cifti_mod(['/home/data/subjects/' subject '/connectome/' functional_sequences{seq} '_corr_smallwall.dconn.nii'],data)
                
                data = [];
                
                
                surface_parcellation_singlesub(subject,['/home/data/subjects/' subject '/connectome/' functional_sequences{seq} '_corr_smallwall.dconn.nii'],parcellation_subsample,0,parcellationfolder)
                
                delete(['/home/data/subjects/' subject '/connectome/' functional_sequences{seq} '_corr_smallwall.dconn.nii'])
                
                ciftifile_normwall = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                parcel_creator_cifti(['corrofcorr_allgrad_LR_subcort_smooth' num2str(geodesic_smooth) '_wateredge_avg.dtseries.nii'],[functional_sequences{seq} '_parcels'],parcel_borderthresh,ciftifile_normwall)
                
                
            end
        end
        
         %% Homogeneity testing
        
        if homogeneity_testing
            mkdir(homogeneityfolder)
            cd(homogeneityfolder)
            dlmwrite(logfile2,'testing homogeneity','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': testing homogeneity against rotated null']);
                
                if ~exist(['/home/data/subjects/' subject '/connectome/' functional_sequences{seq} '_corr.dconn.nii'],'file')
                    ciftifile = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                    data = ft_read_cifti_mod(ciftifile);
                    
                    corr = paircorr_mod(data.data');
                    corr = FisherTransform(corr);
                    corr(isnan(corr)) = 0;
                    
                    data.dimord = 'pos_pos';
                    data.data = corr;
                    dconn = data;
                    clear corr data
                    
                else
                    dconn = ['/home/data/subjects/' subject '/connectome/' functional_sequences{seq} '_corr.dconn.nii'];
                end
                
                outstem = [functional_sequences{seq} '_cov_corr'];
                make_cov_corr(dconn,outstem)
                
                clear dconn
                
                 parcelfile = [parcellationfolder '/' functional_sequences{seq} '_parcels_edgethresh_' num2str(parcel_borderthresh) '.dtseries.nii'];
                
                parcellations_totest = {parcelfile};
                %parcellations_totest = {parcelfile, groupparcelfile};
                
                cov_corr_files = {[outstem '_L.dconn.nii'],[outstem '_R.dconn.nii']};
                
                outdir = homogeneityfolder;
                baddataname = '/home/data/scripts/Resources/Baddata_bigcluster_LR.dtseries.nii';
                
                Batch_homogeneity_testing_cifti(parcellations_totest,cov_corr_files,outdir,[],homogeneity_testing_iterations,baddataname)
                %Batch_homogeneity_testing_cifti_realonly(parcellations_totest,cov_corr_files,outdir,[],homogeneity_testing_iterations,baddataname)
                
                delete(cov_corr_files{1})
                delete(cov_corr_files{2})
                
                parcels = ft_read_cifti_mod(parcelfile);
                homogeneity = ft_read_cifti_mod([homogeneityfolder '/PCA_eigval_per_first_' functional_sequences{seq} '_parcels_edgethresh_' num2str(parcel_borderthresh) '.dtseries.nii']);
                parcels.data(homogeneity.data<homogeneitythresh) = 0;
                ft_write_cifti_mod([parcelfile(1:end-13) '_homogenous.dtseries.nii'],parcels)
                
            end
            
            
        end
                
                
        %% Parcel Infomap and Hub Detection
        
        if parcel_infomap_and_hub_detection
            
           %ign = rmdir(parcel_infomapfolder,'s');
            mkdir(parcel_infomapfolder)
            cd(parcel_infomapfolder)
            
            dlmwrite(logfile2,'running parcel infomap and hub detection','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': running parcel infomap and hub detection']);
                
                parcelfile = [parcellationfolder '/' functional_sequences{seq} '_parcels_edgethresh_' num2str(parcel_borderthresh) '_homogenous.dtseries.nii'];
                %parcelfile = '/home/data/atlases/Group_parcellation/Parcels_LR.dtseries.nii';
%                 cifti_timeseries = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
%                 %cifti_distance_matrix = [ciftifolder '/distances/normalwall_distmat_uint8.mat'];
%                 
%                 if subnum==1
%                     cifti_distance_matrix = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_standardsubcort_uint8.mat';
%                     cifti_distance_matrix = smartload(cifti_distance_matrix);
%                 end
%                 
%                 parcel_infomap_singlesub(parcelfile,cifti_timeseries,cifti_distance_matrix,parcel_infomapfolder,parcel_infomapthresholds,parcel_infomap_xdist,parcel_infomap_minsize)
%                 
%                 consensus_maker_knowncolors([parcel_infomapfolder '/rawassn_minsize' num2str(parcel_infomap_minsize) '_regularized.dscalar.nii'])
% %                consensus_maker_knowncolors_yeo([parcel_infomapfolder '/rawassn_minsize' num2str(parcel_infomap_minsize) '_regularized.dscalar.nii'])

                Hub_detection_summedweightedPC(parcelfile,'corrmat.mat','parcel_distances.mat',[parcel_infomapfolder '/rawassn_minsize' num2str(parcel_infomap_minsize) '.txt'],[parcel_infomapfolder '/rawassn_minsize' num2str(parcel_infomap_minsize) '_regularized_recolored.dscalar.nii'],parcel_infomapthresholds,parcel_infomap_xdist);
                
                
            end
        end
        
        
        
        %% Infomap of Group-defined parcels and Hub Detection
        
        if groupparcel_infomap_and_hub_detection
            
           %ign = rmdir(parcel_infomapfolder,'s');
            mkdir(groupparcel_infomapfolder)
            cd(groupparcel_infomapfolder)
            
            dlmwrite(logfile2,'running group parcel infomap and hub detection','-append','delimiter','')
            
            for seq = 1:length(functional_sequences)
                
                
                disp(['Subject ' subject ', ' functional_sequences{seq} ': running group parcel infomap and hub detection']);
                
                parcelfile = '/home/data/atlases/Group_parcellation/Parcels_LR.dtseries.nii';
                parcelcommunitiesfile = '/home/data/atlases/Group_parcellation/Parcel_Communities.dtseries.nii';
                cifti_timeseries = [ciftifolder '/cifti_timeseries_normalwall/' functional_sequences{seq} '_LR_surf_subcort_333_32k_fsLR_smooth' num2str(geodesic_smooth) '.dtseries.nii'];
                %cifti_distance_matrix = [ciftifolder '/distances/normalwall_distmat_uint8.mat'];
                
%                 if subnum==1
%                     cifti_distance_matrix = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_standardsubcort_uint8.mat';
%                     cifti_distance_matrix = smartload(cifti_distance_matrix);
%                 end
%                 
%                 parcel_infomap_singlesub(parcelfile,cifti_timeseries,cifti_distance_matrix,groupparcel_infomapfolder,parcel_infomapthresholds,parcel_infomap_xdist,parcel_infomap_minsize)
                
                Hub_detection_summedPC(parcelfile,'corrmat.mat','parcel_distances.mat',[groupparcel_infomapfolder '/rawassn_minsize' num2str(parcel_infomap_minsize) '.txt'],parcelcommunitiesfile,parcel_infomapthresholds,parcel_infomap_xdist);
                
                
            end
        end
        
%         %% Hub Detection
%         
%         if hub_detection
%             cd(parcel_infomapfolder)
%             dlmwrite(logfile2,'running hub detection','-append','delimiter','')
%             for seq = 1:length(functional_sequences)
%                 
%                 disp(['Subject ' subject ', ' functional_sequences{seq} ': running hub detection']);
%                 
%                 %parcelfile = [parcellationfolder '/' functional_sequences{seq} '_parcels_edgethresh_' num2str(parcel_borderthresh) '_homogenous.dtseries.nii'];
%                 parcelfile = '/home/data/atlases/Group_parcellation/Parcels_LR.dtseries.nii';
%                 Hub_detection_summedPC(parcelfile,'corrmat.mat','parcel_distances.mat',[parcel_infomapfolder '/rawassn_minsize' num2str(parcel_infomap_minsize) '.txt'],[parcel_infomapfolder '/rawassn_minsize' num2str(parcel_infomap_minsize) '_regularized_recolored.dscalar.nii'],parcel_infomapthresholds,parcel_infomap_xdist);
%                 
%                 
%             end
%         end
        
        
                
            
            
        %%
        
        
        disp(['Subject ' subject ' COMPLETED analysis without error!'])
        dlmwrite(logfile2,'all analysis complete.','-append','delimiter','')
        

    
    
end


