function ComputePC_SingleSubject_medusa(subid,distances,xdistance,tmask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
costs = struct('start',80,'iters',20);
timecoursedata = ft_read_cifti_mod(['/share/cjv2/data/' subid '/fs_LR_output_directory/' subid '/' subid '.Day1.dtseries.nii']);
tmask = load(tmask);
if ischar(distances)
    distances = smartload(distances);
    distances = uint8(distances);
end
%% Load Parcels
parcels_orig = ft_read_cifti_mod(['/share/cjv2/data/' subid '/parcels/Day1/' subid '_parcels_Day1_edgethresh_0.5.dtseries.nii']);
parcels = parcels_orig.data;
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];
parcel_centroids = zeros(length(parcelIDs),1);
parcel_timecourses = zeros(size(timecoursedata.data,2),length(parcelIDs));
for p = 1:length(parcelIDs)
    parcelinds = find(parcels==parcelIDs(p));
    [~,centroidind] = min(sum(distances(parcelinds,parcelinds),2));
    parcel_centroids(p) = parcelinds(centroidind);
    parcel_timecourses(:,p) = mean(timecoursedata.data(parcelinds,:),1);
    parcel_timecourses(:,p) = mean(timecoursedata.data(parcelinds,:),1);
end
          
%% Define Ci derived from the template matching procedures
networks = ft_read_cifti_mod(['/share/cjv2/data/' subid '/parcels/Day1/Networks_Day1_templatematch_dice_kden0.05_parcels.dtseries.nii']);          
%% Create distance exclusions. Default is 20. 
threshdistance = distances > xdistance;
clear distances
threshdistance = threshdistance(parcel_centroids,parcel_centroids);
%% Censor bad time points       
parcel_timecourses(tmask==0,:)=[];

    %% loop through costs
         for cost_thresh_i = 1:costs.iters;
            %% Create correlation matrix using parcel timecourses 
            corrmatrix = paircorr_mod(parcel_timecourses);
            corrmatrix(isnan(corrmatrix)) = 0;
            corrmatrix = FisherTransform(corrmatrix);
            threshmatrix = corrmatrix;
            % Local connections and connections less than the specified cost are set to zero
            threshmatrix(threshdistance==false) = 0;
            threshmatrix(corrmatrix > prctile(icatb_mat2vec(corrmatrix),((costs.start-1)+cost_thresh_i))) = 1;
            corrmatrix(threshmatrix ~= 1) = 0;

            n = networks.data;
            p = parcels_orig.data;
            IDs = unique(p);
 
            %% Set up for Graph Measures based on Template Ci
            for i = 1:length(IDs)
            idx = find(p==IDs(i));
            NetLabel(i) = n(idx(1));
            end
            ZeroPCval = IDs(NetLabel==0);
            NetLabel(IDs<1) = [];
            IDs(IDs<1) = [];
            toremove = find(NetLabel==0);
            IDs(toremove) = []; NetLabel(toremove) = []; corrmatrix(toremove,:) = []; corrmatrix(:,toremove) = [];
            
            %% Correct for networks that are "skipped" during template matching 
            unique_nets = unique(NetLabel);
            for i = 1:length(unique_nets)
                NetLabel(NetLabel==unique_nets(i)) = i;
            end  
            
            %% Replace non-zero values with 1's 
            corrmatrix = icatb_mat2vec(corrmatrix);
            corrmatrix(corrmatrix<prctile(corrmatrix,cost_thresh_i))=0;
            corrmatrix(corrmatrix~=0)=1;
            corrmatrix = icatb_vec2mat(corrmatrix);

            %% Set up for Graph Measures based on Louvain Ci
            NetLabel_louvain = modularity_und(corrmatrix);

            %% Compute all the graph metrics 
    
            PCvals = participation_coef(corrmatrix,NetLabel);
            MDZvals = module_degree_zscore(corrmatrix,NetLabel,0);
            PCvals_louvain = participation_coef(corrmatrix,NetLabel_louvain);
            MDZvals_louvain = module_degree_zscore(corrmatrix,NetLabel_louvain,0);
            
            PCmap = parcels_orig;
            WMDZmap = parcels_orig;
            PCmap_louvain = parcels_orig;
            WMDZmap_louvain = parcels_orig;

            %% Participation Coefficient using the Ci template procedure
            for i = 1:length(IDs)
            PCmap.data(PCmap.data==IDs(i)) = PCvals(i);
            end 
            for i = 1:length(ZeroPCval)
            PCmap.data(PCmap.data==ZeroPCval(i)) = 0;
            end
            %% Participation Coefficient from the Louvain algorithm
            for i = 1:length(IDs)
            PCmap_louvain.data(PCmap_louvain.data==IDs(i)) = PCvals_louvain(i);
            end 
            for i = 1:length(ZeroPCval)
            PCmap_louvain.data(PCmap_louvain.data==ZeroPCval(i)) = 0;
            end
            %% Within Module Degree Z-score using the Ci from the template procedure
            for i = 1:length(IDs)
            WMDZmap.data(WMDZmap.data==IDs(i)) = MDZvals(i);
            end
            for i = 1:length(ZeroPCval)
            WMDZmap.data(WMDZmap.data==ZeroPCval(i)) = 0;
            end
             %% Within Module Degree Z-score using the Ci from the Louvain algorithm
            for i = 1:length(IDs)
            WMDZmap_louvain.data(WMDZmap_louvain.data==IDs(i)) = MDZvals_louvain(i);
            end
            for i = 1:length(ZeroPCval)
            WMDZmap_louvain.data(WMDZmap_louvain.data==ZeroPCval(i)) = 0;
            end

            %% Write out graph theory CIFTI Files and Louvain Ci
            ft_write_cifti_mod(['PCmap_Day1_cost' num2str((costs.start-1)+cost_thresh_i)],PCmap);
            ft_write_cifti_mod(['WMDZ_map_Day1_cost' num2str((costs.start-1)+cost_thresh_i)],WMDZmap);
            ft_write_cifti_mod(['PCmap_Louvain_Day1_cost' num2str((costs.start-1)+cost_thresh_i)],PCmap_louvain);
            ft_write_cifti_mod(['WMDZ_Louvain_Day1_cost' num2str((costs.start-1)+cost_thresh_i)],WMDZmap_louvain);
      
         end          
%% Load all the pc maps
pccost_a = ft_read_cifti_mod(['PCmap_Day1_cost80.dtseries.nii']);
pccost_b = ft_read_cifti_mod(['PCmap_Day1_cost81.dtseries.nii']);
pccost1 = ft_read_cifti_mod(['PCmap_Day1_cost80.dtseries.nii']);
pccost2 = ft_read_cifti_mod(['PCmap_Day1_cost81.dtseries.nii']);
pccost3 = ft_read_cifti_mod(['PCmap_Day1_cost82.dtseries.nii']);
pccost4 = ft_read_cifti_mod(['PCmap_Day1_cost83.dtseries.nii']);
pccost5 = ft_read_cifti_mod(['PCmap_Day1_cost84.dtseries.nii']);
pccost6 = ft_read_cifti_mod(['PCmap_Day1_cost85.dtseries.nii']);
pccost7 = ft_read_cifti_mod(['PCmap_Day1_cost86.dtseries.nii']);
pccost8 = ft_read_cifti_mod(['PCmap_Day1_cost87.dtseries.nii']);
pccost9 = ft_read_cifti_mod(['PCmap_Day1_cost88.dtseries.nii']);
pccost10 = ft_read_cifti_mod(['PCmap_Day1_cost89.dtseries.nii']);
pccost11 = ft_read_cifti_mod(['PCmap_Day1_cost90.dtseries.nii']);
pccost12 = ft_read_cifti_mod(['PCmap_Day1_cost91.dtseries.nii']);
pccost13 = ft_read_cifti_mod(['PCmap_Day1_cost92.dtseries.nii']);
pccost14 = ft_read_cifti_mod(['PCmap_Day1_cost93.dtseries.nii']);
pccost15 = ft_read_cifti_mod(['PCmap_Day1_cost94.dtseries.nii']);
pccost16 = ft_read_cifti_mod(['PCmap_Day1_cost95.dtseries.nii']);
pccost17 = ft_read_cifti_mod(['PCmap_Day1_cost96.dtseries.nii']);
pccost18 = ft_read_cifti_mod(['PCmap_Day1_cost97.dtseries.nii']);
pccost19 = ft_read_cifti_mod(['PCmap_Day1_cost98.dtseries.nii']);
pccost20 = ft_read_cifti_mod(['PCmap_Day1_cost99.dtseries.nii']);

%%  Write out the sum and mean maps

pccost_a.data = (pccost1.data + pccost2.data + pccost3.data + pccost4.data + pccost5.data + pccost6.data + pccost7.data + pccost8.data + pccost9.data + ...
    pccost10.data + pccost11.data + pccost12.data + pccost13.data + pccost14.data + pccost15.data + pccost16.data + pccost17.data + pccost18.data + pccost19.data + pccost20.data);

SumPC = pccost_a;

ft_write_cifti_mod(['PCmap_Day1_sumcost'],SumPC);

pccost_b.data = ((SumPC.data)/20);

MeanPC = pccost_b;

ft_write_cifti_mod(['PCmap_Day1_meancost'],MeanPC);

%% Load all the Within-module z-score maps
wmdzcost_a = ft_read_cifti_mod(['WMDZ_map_Day1_cost80.dtseries.nii']);
wmdzcost_b = ft_read_cifti_mod(['WMDZ_map_Day1_cost81.dtseries.nii']);
wmdzcost1 = ft_read_cifti_mod(['WMDZ_map_Day1_cost80.dtseries.nii']);
wmdzcost2 = ft_read_cifti_mod(['WMDZ_map_Day1_cost81.dtseries.nii']);
wmdzcost3 = ft_read_cifti_mod(['WMDZ_map_Day1_cost82.dtseries.nii']);
wmdzcost4 = ft_read_cifti_mod(['WMDZ_map_Day1_cost83.dtseries.nii']);
wmdzcost5 = ft_read_cifti_mod(['WMDZ_map_Day1_cost84.dtseries.nii']);
wmdzcost6 = ft_read_cifti_mod(['WMDZ_map_Day1_cost85.dtseries.nii']);
wmdzcost7 = ft_read_cifti_mod(['WMDZ_map_Day1_cost86.dtseries.nii']);
wmdzcost8 = ft_read_cifti_mod(['WMDZ_map_Day1_cost87.dtseries.nii']);
wmdzcost9 = ft_read_cifti_mod(['WMDZ_map_Day1_cost88.dtseries.nii']);
wmdzcost10 = ft_read_cifti_mod(['WMDZ_map_Day1_cost89.dtseries.nii']);
wmdzcost11 = ft_read_cifti_mod(['WMDZ_map_Day1_cost90.dtseries.nii']);
wmdzcost12 = ft_read_cifti_mod(['WMDZ_map_Day1_cost91.dtseries.nii']);
wmdzcost13 = ft_read_cifti_mod(['WMDZ_map_Day1_cost92.dtseries.nii']);
wmdzcost14 = ft_read_cifti_mod(['WMDZ_map_Day1_cost93.dtseries.nii']);
wmdzcost15 = ft_read_cifti_mod(['WMDZ_map_Day1_cost94.dtseries.nii']);
wmdzcost16 = ft_read_cifti_mod(['WMDZ_map_Day1_cost95.dtseries.nii']);
wmdzcost17 = ft_read_cifti_mod(['WMDZ_map_Day1_cost96.dtseries.nii']);
wmdzcost18 = ft_read_cifti_mod(['WMDZ_map_Day1_cost97.dtseries.nii']);
wmdzcost19 = ft_read_cifti_mod(['WMDZ_map_Day1_cost98.dtseries.nii']);
wmdzcost20 = ft_read_cifti_mod(['WMDZ_map_Day1_cost99.dtseries.nii']);

%% Write out mean and summed maps

wmdzcost_a.data = (wmdzcost1.data + wmdzcost2.data + wmdzcost3.data + wmdzcost4.data + wmdzcost5.data + wmdzcost6.data + wmdzcost7.data + wmdzcost8.data + wmdzcost9.data + ...
    wmdzcost10.data + wmdzcost11.data + wmdzcost12.data + wmdzcost13.data + wmdzcost14.data + wmdzcost15.data + wmdzcost16.data + wmdzcost17.data + wmdzcost18.data + wmdzcost19.data + wmdzcost20.data);

SumWMDZ = wmdzcost_a;

ft_write_cifti_mod(['WMDZ_map_Day1_sumcost'],SumWMDZ);

wmdzcost_b.data = ((wmdzcost_a.data)/20);

MeanWMDZ = wmdzcost_b;

ft_write_cifti_mod(['WMDZ_map_Day1_meancost'],MeanWMDZ);

end
