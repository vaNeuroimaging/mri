%function Ind_watershed_network_IDs(subjectnum,hem,minimadist,fracmaxh,outputfolder)

% subjectnum = 1;
hem = 'L';
minimadist = 3;
fracmaxh = 1;
outputfolder = '/data/cn4/evan/RestingState/Ind_variability/';
maxdistaway = 50;


medialmaskdata = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']);
medialmaskdata = medialmaskdata.cdata;
corticalindices = find(medialmaskdata==0);
medialindices = find(medialmaskdata);

bufsize=16384;
%caredir = '/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K';
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
neighbors_noNaNs = neighbors;
neighbors_noNaNs(isnan(neighbors_noNaNs)) = medialindices(1);

labelfile = ['/data/cn4/evan/RestingState/Consensus/ConsensusMapvFinal_' hem '.label.gii'];
origlabel = gifti(labelfile);


cohortfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/NEW_nokids_TMASKLIST.txt';
[subjects tmasks] = textread(cohortfile,'%s %s');

disp('Loading average connectivity')
avgconnectivity = gifti(['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/avg_corr_' hem '.func.gii']);
avgconnectivity = FisherTransform(avgconnectivity.cdata);

if strcmp(hem,'L'); hemname = 'LEFT'; elseif strcmp(hem,'R'); hemname = 'RIGHT'; end
roifile = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/' hemname '/cifti_coords.roi'];
[xyz nodenames] = roifilereader(roifile);

roifile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT/cifti_coords.roi';
dmat = euclidean_distance(roifile);
dmat=(dmat>maxdistaway);

numgoodnetworks = length(find(origlabel.private.label.key>0));
numnetworkcomparisons = numgoodnetworks*(numgoodnetworks-1)/2;
pairwisenetworkborders = zeros(length(medialmaskdata),numnetworkcomparisons);

comparisonnum = 0;
for i = 1:numgoodnetworks
    for j = 1:numgoodnetworks
        if j>i
            comparisonnum = comparisonnum+1;
            networkcomparison{comparisonnum} = [origlabel.private.label.name{j+2} ' vs ' origlabel.private.label.name{i+2}];
        end
    end
end

avgwatershedPCA_firsteigval_vertexwise = zeros(length(corticalindices),1);

for subjectnum = 1:length(subjects)
    
    subject = subjects{subjectnum};
    disp(['Subject ' num2str(subjectnum) ' of ' num2str(length(subjects)) ': ' subject])
    
    edgemap = ['/data/cn4/laumannt/left_hem_edge/' subject '_avg_edge_avg_smooth_' hem '_noalone.func.gii'];
    
    newlabel = origlabel;
    newlabel.cdata = newlabel.cdata .* 0;
    
    slashloc = strfind(edgemap,'/');
    filename = edgemap(slashloc(end)+1:end);
    filenamebase = filename(1:end-9);
    watershedfilename = [filenamebase '.func.gii'];
    
    if ~exist([outputfolder filenamebase 'watershedmerge.func.gii'])
        
        %disp('Making minima')
        %metric_minima(edgemap,minimadist,outputfolder,minimafilename,hem);
        
        disp('Making watersheds')
        watersheddata = watershed_algorithm_merge3(edgemap,outputfolder,filenamebase,hem);
        watersheddata = watersheddata(corticalindices);
        %watershed_algorithm_Evan(edgemap,[outputfolder '/' minimafilename],200,fracmaxh,outputfolder,minimafilename(1:end-9));
        
        disp('Finding watershed centers')
        watershedcenters = find_center_water_parcel([outputfolder filenamebase 'watershedmerge.func.gii'],hem);
        watershedcenters = watershedcenters(corticalindices);
    else
    
    watersheddata = gifti([outputfolder filenamebase 'watershedmerge.func.gii']);
    %watersheddata = gifti('/data/cn4/evan/RestingState/Ind_variability/vc32347_minima1watershedmerge.func.gii');
    watersheddata = watersheddata.cdata(corticalindices);
        
    watershedcenters = gifti([outputfolder filenamebase 'watershedmerge_centers.func.gii']);
    %watershedcenters = gifti('/data/cn4/evan/RestingState/Ind_variability/vc32347_minima1watershedmerge_centers.func.gii');
    watershedcenters = watershedcenters.cdata(corticalindices);
    
    end
    
    %delete([outputfolder '/' minimafilename])
    
    subjectdata = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/' subject '_BOLD_L_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii'];
    
    copyfile(subjectdata,'/data/cn4/evan/Temp/Temp.dtseries.nii');
    
    subjectdata = cifti_read('/data/cn4/evan/Temp/Temp.dtseries.nii');
    
    tmask = load(tmasks{subjectnum});
    
    subjectdata = subjectdata(:,logical(tmask));
    
    watersheds = unique(watersheddata);
    watersheds(watersheds==0) = [];
    
    corticallabels = origlabel.cdata(corticalindices);
    
    outputcorticalID = zeros(length(corticalindices),1);
    outputcorticalcorrelation = zeros(length(corticalindices),1);
    outputcorticaldist = zeros(length(corticalindices),1);
    
    disp('Calculating connectivity matrix')
    subcorrelmat = FisherTransform(paircorr_mod(subjectdata'));
    subcorrelmat(isnan(subcorrelmat)) = 0;
    for watershednum = 1:length(watersheds)
        
        watershed = watersheds(watershednum);
        watershedindices = find(watersheddata==watershed);
        
        watershedconnectivity(:,watershednum) = mean(subcorrelmat(:,watershedindices),2);
    end
    allwatershed_correlvsvertices = paircorr_mod(watershedconnectivity,avgconnectivity);
    
    disp('Calculating watershed matches')
    for watershednum = 1:length(watersheds)
        
        string{watershednum} = ['   Watershed ' num2str(watershednum) ' of ' num2str(length(watersheds))];
        if watershednum==1; fprintf('%s',string{watershednum}); else fprintf([repmat('\b',1,length(string{watershednum-1})) '%s'],string{watershednum}); end
        
        watershed = watersheds(watershednum);
        
        watershedcenter = find(watershedcenters==watershed);
        distancevector = dmat(watershedcenter,:);
        
        watershedindices = find(watersheddata==watershed);
        %watershedtimecourse = mean(subjectdata(watershedindices,:),1);
        %watershedtimecourse(isnan(watershedtimecourse)) = 0;
        
        %watershedcorrelmaps = FisherTransform(paircorr_mod(subjectdata(watershedindices,:)',subjectdata'));
        %watershedcorrelmaps(isnan(watershedcorrelmaps)) = 0;
        watershedcorrelmaps = subcorrelmat(watershedindices,:);
        
        if size(watershedcorrelmaps,1) > 2
            
            [ign ign2 eigvals_per] = PCA_reduction(watershedcorrelmaps','comps',2);
            eigval_per_first = eigvals_per(1);
            
        else
            eigval_per_first = 0;
        end
        
        avgwatershedPCA_firsteigval_vertexwise(watershedindices) = avgwatershedPCA_firsteigval_vertexwise(watershedindices) + eigval_per_first;
        
        %watershedcorrelmap = mean(watershedcorrelmaps,1);
        
        %watershed_correlvsvertices = paircorr_mod(watershedcorrelmap',avgconnectivity);
        watershed_correlvsvertices = allwatershed_correlvsvertices(watershednum,:);
        watershed_correlvsvertices(logical(distancevector)) = 0;
        watershed_correlvsvertices = watershed_correlvsvertices(1:length(corticalindices));
        [maxval maxi] = max(watershed_correlvsvertices);
        watershedID = corticallabels(maxi);
        
        watershedxyz = single(xyz(watershedindices,:));
        targetvertxyz = single(xyz(maxi,:));
        distances = sqrt((watershedxyz(:,1) - targetvertxyz(1)).^2 + (watershedxyz(:,2) - targetvertxyz(2)).^2 + (watershedxyz(:,3) - targetvertxyz(3)).^2);
        meandist = mean(distances);
        
        outputcorticalID(watershedindices) = watershedID;
        outputcorticalcorrelation(watershedindices) = maxval;
        outputcorticaldist(watershedindices) = meandist;
        
    end
    
    disp(' ')
    
    newlabel.cdata(corticalindices) = int32(outputcorticalID);
    save(newlabel,[outputfolder '/' subject '_watershed_vertexmatch_networkIDs_merge.label.gii']);
    
    borderindices = intersect(find(newlabel.cdata==0),corticalindices);
    
    disp('   Calculating network pairwise borders')
    
    comparisonnum = 0;
        for network1 = 1:numgoodnetworks
            for network2 = 1:numgoodnetworks
                if network2 > network1
                    comparisonnum = comparisonnum +1;
                    
%                     string{comparisonnum} = networkcomparison{comparisonnum};
%                     if comparisonnum==1; fprintf('%s',string{comparisonnum}); else fprintf([repmat('\b',1,length(string{comparisonnum-1})) '%s'],string{comparisonnum}); end
                    
                    bordersbewteeenthesenetworks = logical((newlabel.cdata==0) .* (medialmaskdata==0) ...
                        .* ((newlabel.cdata(neighbors_noNaNs(:,2))==network1) + (newlabel.cdata(neighbors_noNaNs(:,3))==network1) + (newlabel.cdata(neighbors_noNaNs(:,4))==network1) + (newlabel.cdata(neighbors_noNaNs(:,5))==network1) + (newlabel.cdata(neighbors_noNaNs(:,6))==network1) + (newlabel.cdata(neighbors_noNaNs(:,7))==network1)) ...
                        .* ((newlabel.cdata(neighbors_noNaNs(:,2))==network2) + (newlabel.cdata(neighbors_noNaNs(:,3))==network2) + (newlabel.cdata(neighbors_noNaNs(:,4))==network2) + (newlabel.cdata(neighbors_noNaNs(:,5))==network2) + (newlabel.cdata(neighbors_noNaNs(:,6))==network2) + (newlabel.cdata(neighbors_noNaNs(:,7))==network2)));
                    pairwisenetworkborders(bordersbewteeenthesenetworks,comparisonnum) = pairwisenetworkborders(bordersbewteeenthesenetworks,comparisonnum)+1;
                
                end
            end
        end
        
    
    
    
%     for vert = borderindices';
%         vertneighs = neighbors(vert,2:7);
%         vertneighs(isnan(vertneighs)) = [];
%         
%         comparisonnum = 0;
%         for network1 = 1:numgoodnetworks
%             for network2 = 1:numgoodnetworks
%                 if network2 > network1
%                     comparisonnum = comparisonnum +1;
%                     
%                     if any(newlabel.cdata(vertneighs)==network1) && any(newlabel.cdata(vertneighs)==network2)
%                         pairwisenetworkborders(vert,comparisonnum) = pairwisenetworkborders(vert,comparisonnum)+1;
%                     end
%                 end
%             end
%         end
%     end
    
    
    % output = zeros(length(newlabel.cdata),1);
    % output(corticalindices) = outputcorticalID;
    % save(gifti(single(output)),[outputfolder '/' subject '_watershed_vertexmatch_networkIDs.func.gii']);
    
    
%     output = zeros(length(newlabel.cdata),1);
%     output(corticalindices) = outputcorticalcorrelation;
%     save(gifti(single(output)),[outputfolder '/' subject '_watershed_vertexmatch_correlstrength.func.gii']);
%     
%     output = zeros(length(newlabel.cdata),1);
%     output(corticalindices) = outputcorticaldist;
%     save(gifti(single(output)),[outputfolder '/' subject '_watershed_vertexmatch_dist.func.gii']);
    
end

pairwisenetworkborders_output = gifti(single(pairwisenetworkborders));
attributesstuff = pairwisenetworkborders_output.private.data{1}.attributes;
attributesstuff.Dim = [size(pairwisenetworkborders,1) 1];
for i = 1:size(pairwisenetworkborders,2)
    pairwisenetworkborders_output.private.data{i}.attributes = attributesstuff;
    pairwisenetworkborders_output.private.data{i}.data = single(pairwisenetworkborders(:,i));
    pairwisenetworkborders_output.private.data{i}.metadata(1).name = 'Name';
    pairwisenetworkborders_output.private.data{i}.metadata(1).value = networkcomparison{i};
    pairwisenetworkborders_output.private.data{i}.space = [];
end

save(pairwisenetworkborders_output,[outputfolder 'Pairwise_network_border_counts.func.gii']);


avgwatershedPCA_firsteigval_vertexwise_output = zeros(length(newlabel.cdata),1);
avgwatershedPCA_firsteigval_vertexwise_output(corticalindices) = avgwatershedPCA_firsteigval_vertexwise ./ length(subjects);

save(gifti(single(avgwatershedPCA_firsteigval_vertexwise_output)),[outputfolder 'watershedPCA_avgfirsteigval.func.gii'])