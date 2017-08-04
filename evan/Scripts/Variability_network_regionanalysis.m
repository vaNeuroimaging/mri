
clustersizefileL = '/data/cn4/evan/RestingState/Ind_variability/120/Variability_L_consensus_dice_clustersize.func.gii';
clustersizefileR = '/data/cn4/evan/RestingState/Ind_variability/120/Variability_R_consensus_dice_clustersize.func.gii';

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;


%gifti_to_cifti('Variabile_regions_L_consensus_dice_distance10.func.gii','Variabile_regions_R_consensus_dice_distance10.func.gii','Variabile_regions_LR_consensus_dice_distance10')

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');
for i = 1:length(tmasks)
tmask = load(tmasks{i});
ntimepoints(i) = nnz(tmask);
end

all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance10.dtseries.nii');

  networkIDtotest = 7;
  all_variable_regions(all_variable_regions~=networkIDtotest) = 0;

variable_regions = metric_cluster_cifti(all_variable_regions,.5,100,0);
cifti_write_wHDR(variable_regions,[],'Variabile_regions_LR_consensus_dice_distance10_separated_VA.dtseries.nii')
% 
% 
% cifti_to_gifti('Variabile_regions_LR_consensus_dice_distance10_separated.dtseries.nii');
% hems = {'L','R'};
% for hemnum = 1:length(hems)
%     hem = hems{hemnum};
%     system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -metric-resample Variabile_regions_LR_consensus_dice_distance10_separated_' hem '.func.gii /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.' hem '.sphere.164k_fs_LR.surf.gii BARYCENTRIC  Variabile_regions_LR_consensus_dice_distance10_separated_' hem '_164.func.gii -largest'])
%     variable_regions_164 = gifti(['Variabile_regions_LR_consensus_dice_distance10_separated_' hem '_164.func.gii']); variable_regions_164 = variable_regions_164.cdata;
%     fullvariabilitymap = gifti(['Variability_' hem '_distantonly_164.func.gii']); fullvariabilitymap = fullvariabilitymap.cdata;
%     out = zeros(size(variable_regions_164,1),0);
%     for i = 1:size(variable_regions_164,2)
%         if nnz(variable_regions_164(:,i)) > 0;
%             out(:,end+1) = fullvariabilitymap .* variable_regions_164(:,i);
%         end
%     end
%     save(gifti(single(out)),['Variability_' hem '_distantonly_separated_164.func.gii'])
% end
            

%variable_regions(:,1) = [];



%%

surfaceareas = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');

networkconnection_bysub = cifti_read('Templatematch_dice_bysubject.dtseries.nii');

networkconnections = zeros(size(variable_regions,2),size(networkconnection_bysub,2));
networkconnections_SAincluster = zeros(size(variable_regions,2),size(networkconnection_bysub,2));

for i = 1:size(variable_regions,2)
    
    targetID = mode(variable_regions(logical(variable_regions(:,i)),i));
    
    for s = 1:size(networkconnection_bysub,2)
        
        %networkconnections(i,s) = mode(networkconnection_bysub(logical(variable_regions(:,i)),s));
        
        if sum(surfaceareas(logical((networkconnection_bysub(:,s)==targetID) .* (variable_regions(:,i))))) >= (.25 * sum(surfaceareas(logical(variable_regions(:,i)))))
            %nnz(networkconnection_bysub(logical(variable_regions(:,i)),s)==targetID) >= (.05 * nnz(variable_regions(:,i)))
            %networkconnections(i,s) = targetID;
            networkconnections(i,s) = 1;
        else
            %networkconnections(i,s) = mode(networkconnection_bysub(logical(variable_regions(:,i)),s));
            networkconnections(i,s) = 0;
        end
        
        networkconnections_SAincluster(i,s) = sum(surfaceareas(logical((networkconnection_bysub(:,s)==targetID) .* (variable_regions(:,i)))));
        
    end
    
end



[R,P]=corrcoef(networkconnections_SAincluster');
%disp(R)
%disp(P .* nnz(triu(P,1)))

Pcorrected = P .* nnz(triu(P,1));
sigindices = (Pcorrected < 0.05) .* ceil(triu(P,1));
disp(R .* sigindices)
disp(Pcorrected .* sigindices)

correlatedspots = zeros(size(all_variable_regions,1),0);

[x,y] = find(sigindices);
for i = 1:length(x)
    correlatedspots(:,end+1) = sum(variable_regions(:,[x(i) y(i)]),2);
end

%cifti_write_wHDR(correlatedspots,[],'VariableRegions_Correlated')


%%

percentregion_thresh = .25;

all_variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8_selected.dtseries.nii');
variable_regions = cifti_read('Variabile_regions_LR_consensus_dice_distance8separated_selected.dtseries.nii');
networkconnection_bysub = cifti_read('Templatematch_dice_bysubject_kden0.05.dtseries.nii');
mostcommon = cifti_read('Variability_LR_consensus_dice_mostcommon.dtseries.nii');

surfdatafile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
[subjects ciftifiles] = textread(surfdatafile,'%s %s');

tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

ncortverts = 59412;

%main = zeros(ncortverts,size(variable_regions,2));
%alternate = zeros(ncortverts,size(variable_regions,2));
alternatecount = zeros(1,size(variable_regions,2));

for s = 1:length(subjects)
    disp(subjects{s})
    tmask = load(tmasks{s});
    subdata = cifti_read(ciftifiles{s});
    subdata = subdata(:,logical(tmask));
    
    for r = 1:size(variable_regions,2)
        
        if s==1
            alternate{r} = zeros(ncortverts,0);
            main{r} = zeros(ncortverts,0);
        end
        
        mainID = mode(mostcommon(logical(variable_regions(:,r)))); 
        altID = mean(all_variable_regions(logical(variable_regions(:,r))));
        if (nnz(networkconnection_bysub(logical(variable_regions(:,r)),s)==altID) / nnz(variable_regions(:,r))) >= percentregion_thresh
            alternatecount(r) = alternatecount(r)+1;
            
            indices = logical(variable_regions(:,r) .* (networkconnection_bysub(:,s)==altID));
            
            corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
            corrpattern(isnan(corrpattern)) = 0;
            
            %alternate(:,r) = alternate(:,r) + corrpattern(1:ncortverts);
            alternate{r}(:,end+1) = corrpattern(1:ncortverts);
            
        else
            
            indices = logical(variable_regions(:,r) .* (networkconnection_bysub(:,s)==mainID));
            
            corrpattern = paircorr_mod(subdata',mean(subdata(indices,:),1)');
            corrpattern(isnan(corrpattern)) = 0;
            
            %main(:,r) = main(:,r) + corrpattern(1:ncortverts);
            main{r}(:,end+1) = corrpattern(1:ncortverts);
        end
    end
end

mainout = zeros(size(mostcommon,1),size(variable_regions,2));
altout = zeros(size(mostcommon,1),size(variable_regions,2));
tout = zeros(size(mostcommon,1),size(variable_regions,2));
tthresh = tinv(1-(.05 / ncortverts / size(variable_regions,2)),length(subjects));
for r = 1:size(variable_regions,2)
    
    %     alternate(:,r) = alternate(:,r) ./ alternatecount(r);
    %     main(:,r) = main(:,r) ./ (length(subjects) - alternatecount(r));
    
    mainout(1:ncortverts,r) = mean(main{r},2);
    altout(1:ncortverts,r) = mean(alternate{r},2);
    
    [H,P,CI,STATS] = ttest2(main{r}',alternate{r}');
    tout(1:ncortverts,r) = STATS.tstat;
    
end

disp(alternatecount)

cifti_write_wHDR(mainout,[],'MeanConnectivity_mainID_byregion')
cifti_write_wHDR(altout,[],'MeanConnectivity_alternateID_byregion')
cifti_write_wHDR(tout,[],['MeanConnectivity_T_MainVsAlternateID_byregion_correctedT' num2str(tthresh)])
cifti_write_wHDR((mainout - altout),[],'MeanConnectivity_Diff_MainVsAlternateID_byregion')
        
    
