%%

% articfile = '/data/cn4/evan/RestingState/Ind_variability/articpoints.txt';
% articdata = load(articfile);
% articcoords = articdata(:,1:3);

articcoords = [6 42 8; 2 13 35];%; -34 16 3; -50 -63 3; -2 58 -8; -2 -50 36];

roifile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/RIGHT/cifti_coords.roi';
[xyz nodenames] = roifilereader(roifile);
xyz = single(xyz);

for i = 1:size(articcoords,1)
    distancevec = sqrt((xyz(:,1)-articcoords(i,1)).^2 + (xyz(:,2)-articcoords(i,2)).^2 + (xyz(:,3)-articcoords(i,3)).^2);
    [minval, index] = min(distancevec);
    disp(['Point ' num2str(i) ' (' num2str(articcoords(i,:)) '): ' num2str(index)])
end



%%
ciftiindex = 27962;

cohortfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/NEW_nokids_TMASKLIST.txt';

medialmaskdata = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');

[subjects tmasks] = textread(cohortfile,'%s %s');

for s = 1:length(subjects)
        
        subject = subjects{s};
        
        string{s} = ['    Subject ' num2str(s) ': ' subject];
        if s==1; fprintf('%s',string{s}); else fprintf([repmat('\b',1,length(string{s-1})) '%s'],string{s}); end
        
        subjectdata = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/' subject '_BOLD_L_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii'];
        
        copyfile(subjectdata,'/data/cn4/evan/Temp/Temp.dtseries.nii');
        
        subjectdata = cifti_read('/data/cn4/evan/Temp/Temp.dtseries.nii');
        
        tmask = load(tmasks{s});
        
        subjectdata = subjectdata(:,logical(tmask))';
        
        correlmaps(s,:) = paircorr_mod(subjectdata(:,ciftiindex),subjectdata)';

        names{1,s} = num2str(s);
        
end

%correlmaps(isnan(correlmaps)) = 0;

nanvals = single(isnan(correlmaps));
nancols = logical(sum(nanvals,1));

origcorrelmaps = correlmaps;

correlmaps(:,nancols) = [];

Y = pdist(correlmaps,'correlation');

clustering = linkage(Y, 'average');

[H,T,perm] = dendrogram(clustering, 0, 'orientation','left','labels', names, 'colorthreshold', .9);         % 'dendrogram' creates the tree
orient landscape;   
clusters = cluster(clustering, 'Cutoff', 0.5, 'Criterion', 'distance');										% 'cluster' reorders the regions as they 																												  appear on the dendrogram
       
cophenetic_r = cophenet(clustering, Y)
%%
%clusters = {[1:46],[48:53],[54:56],[57:63],[65:105],[106:120]};
%clusters = {[1:53],[54:63],[65:105],[106:120]};
%clusters = {[1:64],[65:120]};
clusters = {[1:86],[93:101],[102:111],[113:117]};

output = zeros(size(medialmaskdata.cdata,1),length(clusters));

for cluster = 1:length(clusters)
    
    meanclustermap = mean(origcorrelmaps(perm(clusters{cluster}),:),1)';
    
    clustercortex = meanclustermap(1:length(find(medialmaskdata.cdata==0)));
    
    output(~medialmaskdata.cdata,cluster) = clustercortex;
    
    
end

save(gifti(single(output)),['/data/cn4/evan/RestingState/Ind_variability/mypMFG_Clustering_' num2str(length(clusters)) 'clusters_MeanConnectivity.func.gii'])
