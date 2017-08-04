eigval_perc_notfirsttwooutputfilename = '/data/cn4/evan/RestingState/Ind_variability/eigval_perc_notfirsttwo_vertexwise_L.func.gii';
numcomponents_30percoutputfilename = '/data/cn4/evan/RestingState/Ind_variability/numcomponents_30perc_vertexwise_L.func.gii';

cohortfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/NEW_nokids_TMASKLIST.txt';

xdistance = 20;

medialmaskdata = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');

divisions = 20;

roifile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT/cifti_coords.roi';
dmat = euclidean_distance(roifile);
dmat=(dmat<xdistance);

nodesperdivision = ceil(size(dmat,1) / divisions);

[subjects tmasks] = textread(cohortfile,'%s %s');
warning off

eigval_perc_notfirsttwo = zeros(size(dmat,1),1);
numcomponents_30perc = zeros(size(dmat,1),1);
%first4components = zeros(size(allsubjectdata{1},1),size(allsubjectdata{1},1),4);

for divisionnum = 1:divisions
    numallsubjects = 0;
    %get the vertex indices for this division
    if divisionnum==divisions
        indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : size(dmat,1);
    else
        indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : nodesperdivision*divisionnum;
    end

    
    disp(['Division ' num2str(divisionnum)])
    
    divisionfilename{divisionnum} = ['/data/cn4/evan/RestingState/Ind_variability/Division' num2str(divisionnum) '.mat'];
    
    if ~exist(divisionfilename{divisionnum})
    
    for s = 1:length(subjects)
        
        subject = subjects{s};
        
        string{s} = ['Subject ' num2str(s) ': ' subject];
        if s==1; fprintf('%s',string{s}); else fprintf([repmat('\b',1,length(string{s-1})) '%s'],string{s}); end
        
        subjectdata = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/' subject '_BOLD_L_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii'];
        
        copyfile(subjectdata,'/data/cn4/evan/Temp/Temp2.dtseries.nii');
        
        subjectdata = cifti_read('/data/cn4/evan/Temp/Temp2.dtseries.nii');
        
        tmask = load(tmasks{s});
        
        subjectcorrel = paircorr_mod(single(subjectdata(indices{divisionnum},:)'),single(subjectdata'));
        
        subjectcorrelinindices = FisherTransform(subjectcorrel);
        
        valuestozero = (isnan(subjectcorrelinindices)) + dmat(indices{divisionnum},:);
        
        subjectcorrelinindices(logical(valuestozero)) = 0;
        
        subjectcorrels(:,:,s) = subjectcorrelinindices;
        
        clear subjectdata subjectcorrelinindices
        
    end
    
    %save(divisionfilename{divisionnum},'subjectcorrels','-v7.3')
    
    disp(' ')
    
    else
        
        load(divisionfilename{divisionnum})
    
    end
    
    
    
    for grord = 1:length(indices{divisionnum})
        
        string{grord} = ['GrayOrdinate ' num2str(grord) ' of ' num2str(length(indices{divisionnum}))];
        if grord==1; fprintf('%s',string{grord}); else fprintf([repmat('\b',1,length(string{grord-1})) '%s'],string{grord}); end
        
        [Y(:,1,:) V_use eigvals_per] = PCA_reduction(squeeze(subjectcorrels(grord,:,:)),'comps',8);
        
        eigval_perc_notfirsttwo(indices{divisionnum}(grord)) = 100-(eigvals_per(1)+eigvals_per(2));
        
        totalperc = 0;
        count = 0;
        while totalperc < 30
            count = count+1;
            totalperc = totalperc + eigvals_per(count);
        end
        numcomponents_30perc(indices{divisionnum}(grord)) = count;
        
        %first4components(:,grord,:) = Y(:,1,1:4);
        
        clear Y
        
    end
    
    disp(' ')
    
    clear subjectcorrels
    
end


output = zeros(size(medialmaskdata.cdata));

output(~medialmaskdata.cdata) = eigval_perc_notfirsttwo(1:29696);

save(gifti(single(output)),eigval_perc_notfirsttwooutputfilename);

output(~medialmaskdata.cdata) = numcomponents_30perc(1:29696);

save(gifti(single(output)),numcomponents_30percoutputfilename);


