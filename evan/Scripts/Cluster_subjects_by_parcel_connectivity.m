datalist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/FC_Mapping_120/AllC_TMASKLIST.txt';
iscifti = 2;
hem = 'L';

kdenthresholds = [.2 .2];
kdeninterval = .01;
networksizeminimum = 10;

variableparcelIDs = [15024];

outputfolder = ['/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_latOccip/'];

parcels = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/120_L_crossthresh_watershedmerge.func.gii'); parcels = parcels.cdata;

baddata = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); baddata = baddata.cdata;
baddata = baddata<800;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];


% parcelIDs_incrap = unique(parcels .* baddata);
% for badparcelID = parcelIDs_incrap'
%     parcelIDs(parcelIDs == badparcelID) = [];
% end

% subject_correlations_ofparcels = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/Subject_avgcorrelation_to_120_L_crossthresh_watershedmerge.func.gii');
% subject_correlations_ofparcels = subject_correlations_ofparcels.cdata;
 for parcelnum = 1:length(parcelIDs)
%     variableparcelindex(parcelnum) = mean(subject_correlations_ofparcels(parcels==parcelIDs(parcelnum))) < .55;
    variableparcelindex(parcelnum) = any(variableparcelIDs==parcelIDs(parcelnum));
 end


%origparcels = gifti(parcelfilename); origparcels = origparcels.cdata;



%%

if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
else
    parcels = origparcels;
    ncortexLverts = length(parcels);
end


for parcelnum = 1:length(parcelIDs)
    parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
end



[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

allsubcorrelmats = zeros(length(parcelIDs)*nnz(variableparcelindex),length(subjects));

for subnum = 1:length(subjects)
    
    disp(['Subject ' num2str(subnum)])
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{subnum});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    parceltimecourses = zeros(length(parcelIDs),size(subtimecourse,2));
    
    for parcelnum = 1:length(parcelIDs)
        parceltimecourses(parcelnum,:) = mean(subtimecourse(parcelindices{parcelnum},:),1);
    end
    
    parcelcorrelmat = paircorr_mod(parceltimecourses(logical(variableparcelindex),:)',parceltimecourses');
    
    allsubcorrelmats(:,subnum) = reshape(parcelcorrelmat,numel(parcelcorrelmat),1);
    names{subnum} = ['Subject ' num2str(subnum)];
end
 
%%
Y = pdist(allsubcorrelmats','correlation');           % 'pdist' converts the square adjacency matrix to a 
                                         %  1 x n matrix so that the function linkage can construct the tree
                                      
clustering = linkage(Y, 'average');      % 'linkage' computes the data to construct the tree
										 % 'average' refers to the UPGMA algorithm


[H,T,perm] = dendrogram(clustering, 0, 'orientation','left','labels', names, 'colorthreshold', .63);         % 'dendrogram' creates the tree

clusters = cluster(clustering, 'Cutoff', 0.5, 'Criterion', 'distance');										% 'cluster' reorders the regions as they 																												  appear on the dendrogram

orient landscape;                                                         							        % orients the dendrogram to 
																											% either landscape (as shown) or portrait

cophenetic_r = cophenet(clustering, Y)

set(gcf,'Nextplot','new');

%%

sub_similaritymat = paircorr_mod(allsubcorrelmats);

imagesc(sub_similaritymat)
set(gcf,'Nextplot','new');

corrmatname = 'subject_corrmat.mat';
save(corrmatname,'sub_similaritymat');


 roifilename = 'subjects.roi';
quickroifile(length(subjects),roifilename);

prmfilename = 'prmfile.txt';
delete(prmfilename)
fid = fopen([prmfilename],'at'); %open the output file for writing
fclose(fid);

dlmwrite(prmfilename,corrmatname,'-append','delimiter','');
dlmwrite(prmfilename,roifilename,'-append','delimiter','');
dlmwrite(prmfilename,'watershed_bothhem','-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(1)),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdeninterval),'-append','delimiter','');
dlmwrite(prmfilename,num2str(kdenthresholds(2)),'-append','delimiter','');
dlmwrite(prmfilename,outputfolder,'-append','delimiter','');
dlmwrite(prmfilename,num2str(.1),'-append','delimiter','');
dlmwrite(prmfilename,num2str(12),'-append','delimiter','');
dlmwrite(prmfilename,num2str(5),'-append','delimiter','');
dlmwrite(prmfilename,'kden','-append','delimiter','');

%%
graphcluster_Evan(prmfilename,'thr',[],1,0,'infomap');

%% Modify color assignments with miminum network size criterion

filesinoutputfolder = dir(outputfolder);
for i=3:length(filesinoutputfolder);
    folderdatecreated(i-2) = filesinoutputfolder(i).datenum;
end

[maxval maxindex] = max(folderdatecreated);
trueoutputfolder = [outputfolder '/' filesinoutputfolder(maxindex+2).name];

cd(trueoutputfolder)

simple_assigns = modify_clrfile('simplify','rawassn.txt',networksizeminimum);