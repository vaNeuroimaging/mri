datalist = 'AllC_DATALIST.txt';
tmasklist = 'AllC_TMASKLIST.txt';
iscifti = 2;
hem = 'L';

kdenthresholds = [.01 .2];
kdeninterval = .01;
networksizeminimum = 3;

outputfolder = ['/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_SotS_v3/'];

ROIfile = 'PostTempHubL.func.gii';





%origparcels = gifti(parcelfilename); origparcels = origparcels.cdata;



%%
mkdir(outputfolder);
ROI = gifti(ROIfile); ROI = ROI.cdata;

if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    ROI = ROI(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    ROI = ROI(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
else
    ROI = origparcels;
    ncortexLverts = length(ROI);
end

ROIindices = find(ROI) + (strcmp(hem,'R') * ncortexLverts);



[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end


for subnum = 1:length(subjects)
    
    disp(['Subject ' num2str(subnum)])
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' ' outputfolder '/Temp.func.gii']);
    subtimecourse = gifti([outputfolder '/Temp.func.gii']); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{subnum});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    ROItimecourse = mean(subtimecourse(ROIindices,:),1);
    
    ROIconnectivity = paircorr_mod(ROItimecourse',subtimecourse');
    
    allsubcorrelmats(:,subnum) = reshape(ROIconnectivity,numel(ROIconnectivity),1);
    names{subnum} = ['Subject ' subjects{subnum}];
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
% addpath /data/cn4/evan/Scripts/
% graphcluster_Evan(prmfilename,'thr',[],1,0,'infomap');
ngt.graphcluster(prmfilename,'thr',[],1,0,'infomap');

%% Modify color assignments with miminum network size criterion

filesinoutputfolder = dir(outputfolder);
for i=3:length(filesinoutputfolder);
    folderdatecreated(i-2) = filesinoutputfolder(i).datenum;
end

[maxval maxindex] = max(folderdatecreated);
trueoutputfolder = [outputfolder '/' filesinoutputfolder(maxindex+2).name];

cd(trueoutputfolder)

simple_assigns = modify_clrfile('simplify','rawassn.txt',networksizeminimum);