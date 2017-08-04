avgoutputfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/Avg_variability_L.func.gii';
stdoutputfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/Std_variability_L.func.gii';

cohortfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/NEW_nokids_TMASKLIST.txt';

xdistance = 20;

medialmaskdata = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');

divisions = 20;

roifile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT/cifti_coords.roi';
dmat = euclidean_distance(roifile);
dmat=(dmat>=xdistance);

[subjects tmasks] = textread(cohortfile,'%s %s');

correls = zeros(size(dmat,1),(length(subjects)^2 - length(subjects))/2);

nodesperdivision = ceil(size(dmat,1) / divisions);

warning off

disp('Calculating correlation between surface nodes and volume')

%Loop through divisions
for divisionnum = 1:divisions
    numallsubjects = 0;
    %get the vertex indices for this division
    if divisionnum==divisions
        indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : size(dmat,1);
    else
        indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : nodesperdivision*divisionnum;
    end
    
    disp(['Division ' num2str(divisionnum)])
    
    for s = 1:length(subjects)
        
        subject = subjects{s};
        
        string{s} = ['    Subject ' num2str(s) ': ' subject];
        if s==1; fprintf('%s',string{s}); else fprintf([repmat('\b',1,length(string{s-1})) '%s'],string{s}); end
        
        subjectdata = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/' subject '_BOLD_L_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii'];
        
        copyfile(subjectdata,'/data/cn4/evan/Temp/Temp2.dtseries.nii');
        
        subjectdata = cifti_read('/data/cn4/evan/Temp/Temp2.dtseries.nii');
        
        tmask = load(tmasks{s});
        
        subjectdata = subjectdata(:,logical(tmask));
        
        subjectcorrel = paircorr_mod(single(subjectdata'));
        
        subjectcorrelinindices = FisherTransform(subjectcorrel(indices{divisionnum},:));
        
        valuestozero = (isnan(subjectcorrelinindices)) + ~dmat(indices{divisionnum},:);
        
        subjectcorrelinindices(logical(valuestozero)) = 0;
        
        subjectcorrels(:,:,s) = subjectcorrelinindices;
        
        
    end
    
    clear subjectcorrelinindices valuestozero subjectcorrel subjectdata
    
    disp(' ')
    clear string
    
    paircount = 0;
    
    for i = 1:length(subjects)
        for j = 1:length(subjects)
            if i>j
                paircount = paircount+1;
                string{paircount} = ['    Comparison ' num2str(paircount) ' of ' num2str((length(subjects)^2-length(subjects))/2)];
                if paircount==1; fprintf('%s',string{paircount}); else fprintf([repmat('\b',1,length(string{paircount-1})) '%s'],string{paircount}); end
                
                subcomparison = paircorr_mod(squeeze(subjectcorrels(:,:,i))',squeeze(subjectcorrels(:,:,j))');
                allsubcomparisons(:,paircount) = diag(subcomparison);
                
            end
        end
        
    end
    disp(' ')
    
    wholebrainallcomparisons(indices{divisionnum},1:paircount) = allsubcomparisons;
    
    clear subjectcorrel subjectcorrels subjectdata allsubcomparisons
    
end

Zcorrels = FisherTransform(wholebrainallcomparisons);
avgcorrels = mean(Zcorrels,2);
stdcorrels = std(Zcorrels,0,2);

output = zeros(size(medialmaskdata.cdata));

output(~medialmaskdata.cdata) = avgcorrels(1:29696);

save(gifti(single(output)),avgoutputfilename);

output(~medialmaskdata.cdata) = stdcorrels(1:29696);

save(gifti(single(output)),stdoutputfilename);


