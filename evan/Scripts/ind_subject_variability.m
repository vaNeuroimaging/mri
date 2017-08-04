avgoutputfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/Avg_variability_L.func.gii';
stdoutputfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/Std_variability_L.func.gii';

cohortfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/NEW_nokids_TMASKLIST.txt';

xdistance = 20;

medialmaskdata = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');

groupavg = gifti('/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/avg_corr_L.func.gii');
groupavg = groupavg.cdata;
groupavg = FisherTransform(groupavg);

roifile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT/cifti_coords.roi';
dmat = euclidean_distance(roifile);
dmat=(dmat>=xdistance);

[subjects tmasks] = textread(cohortfile,'%s %s');

correls = zeros(size(dmat,1),length(subjects));

for s = 1:length(subjects)
    
    subject = subjects{s};
    
    disp(['Subject ' num2str(s) ': ' subject])
    
    subjectdata = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/' subject '_BOLD_L_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii'];
    
    copyfile(subjectdata,'/data/cn4/evan/Temp/Temp.dtseries.nii');
    
    subjectdata = cifti_read('/data/cn4/evan/Temp/Temp.dtseries.nii');
    
    tmask = load(tmasks{s});
    
    subjectdata = subjectdata(:,logical(tmask));
    
    subjectcorrel = paircorr_mod(single(subjectdata'));
    
    subjectcorrel = FisherTransform(subjectcorrel);
    
    for grord = 1:size(correls,1)
        
        indicestouse = intersect(find(dmat(:,grord)),find(~isnan(subjectcorrel(:,grord))));
        
        if ~isempty(indicestouse)
        
            correls(grord,s) = corr(groupavg(indicestouse,grord),subjectcorrel(indicestouse,grord));
            
        else
            
            correls(grord,s) = 0;
        
        end
        
    end
    
end

Zcorrels = FisherTransform(correls);
avgcorrels = mean(Zcorrels,2);
stdcorrels = std(Zcorrels,0,2);

output = zeros(size(medialmaskdata.cdata));

output(~medialmaskdata.cdata) = avgcorrels(1:29696);

save(gifti(single(output)),avgoutputfilename);

output(~medialmaskdata.cdata) = stdcorrels(1:29696);

save(gifti(single(output)),stdoutputfilename);


