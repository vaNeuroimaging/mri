subjects = textread('/home/data/subjects/processing_list_10min.txt','%s');
subjects(end+1) = {'MAV016'};

ROIs = load_untouch_nii('/home/data/atlases/BB264_MNI_333_all.nii.gz');
ROIs = reshape(ROIs.img,[prod(ROIs.hdr.dime.dim(2:4)) ROIs.hdr.dime.dim(5)]);
values = unique(ROIs); values(values<1) = [];
[x y z] = textread('/home/data/atlases/BB264_MNI_333_all_coords.txt','%n%n%n');

coords_dist = squareform(pdist([x y z]));

coords_dist_short = coords_dist<50;

cross_subject_mean_FD_box = [];
cross_subject_mean_diff_v_rand = [];

FDthresh = 0.05:.03:.5;
    allslopes = zeros(length(subjects),length(FDthresh));

for s = 2:length(subjects)
    subdir = ['/home/data/subjects/' subjects{s} '/preprocessed/'];
    restruns = dir([subdir 'restingstate_*_st_mcf_MNI.nii.gz']);
    motionparams = dir([subdir 'restingstate_*_mcf.par']);
    
    brainmask = load_untouch_nii(['/home/data/subjects/' subjects{s} '/freesurfer/nusmask/aparc+aseg_brainmask_mask_333.nii.gz']);
    brainmask = reshape(brainmask.img,[prod(brainmask.hdr.dime.dim(2:4)) brainmask.hdr.dime.dim(5)]);
    
    tc_concat{s} = zeros(0,length(x));
    FD_concat{s} = zeros(0,1);
    DV_concat{s} = zeros(0,1);
    
    prevstring = [];
    for i = 1:length(restruns)
        string = [subjects{s} ',run ' num2str(i)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        thisrun_params = load([subdir motionparams(i).name]);
        
        thisrun_rot = thisrun_params(:,1:3);
        thisrun_rot_mm = thisrun_rot * 50;
        thisrun_params(:,1:3) = thisrun_rot_mm;
        thisrun_params_delta = [zeros(1,6) ; [thisrun_params(2:end,:) - thisrun_params(1:end-1,:)]];
        thisrun_FD = [sum(abs(thisrun_params_delta),2)];
        
        FD_concat{s} = [FD_concat{s}; thisrun_FD];
        
        
        
         data = load_untouch_nii([subdir restruns(i).name]);
         data = reshape(data.img,[prod(data.hdr.dime.dim(2:4)) data.hdr.dime.dim(5)]);
        gs = mean(data(logical(brainmask),:),1);
        regressors = gs'; %[1:size(data,2)]' ones(size(data,2),1)];
        [residualdata,~,~] = regress_nuisance(data,regressors);
        
        thisrun_DV = [0 abs(diff(gs))];
        DV_concat{s} = [DV_concat{s} thisrun_DV];

        %residualdata = data;
        
        tcs = zeros(size(data,2),length(x));
        
        for roinum = 1:length(values)
            tcs(:,roinum) = mean(residualdata(ROIs==values(roinum),:),1)';
        end
        
        tc_concat{s} = [tc_concat{s}; tcs];
                
        
    end
    disp(' ')
    %% Calculate corrmat from low FD frames
     
    all_corrmat = corrcoef(tc_concat{s});
    
    [FD_concat_sort, ind] = sort(FD_concat{s},'ascend');
    %[DV_concat_sort, ind] = sort(DV_concat{s},'ascend');
    
    mean_diff = [];
    mean_FD_box = [];
    mean_rand_diff = [];
    
    gooddata_induse = find(FD_concat{s} < .1);%ind(1:1000); %true(length(FD_concat{s}),1); %find(DV_concat{s} < 1);
    good_corrmat = corrcoef(tc_concat{s}(gooddata_induse,:));
    
    
    
    %% Calculate differene between correlations from frames at sliding FD to low FD frames
    %binlength = round(min(length(gooddata_induse)/20,500));
    %interval = binlength;
    %binnum = floor((length(FD_concat{s})-binlength)/interval)+1;
    
    binnum = 30;
    binlength = floor(length(FD_concat{s}) / binnum);
    interval = binlength;
    
    mask = ones(length(values));
    mask = triu(mask,1);
    mask = mask & coords_dist_short; %Only look at short distance correlations
    clear mean_FD_box mean_diff mean_rand_diff
    prevstring = [];
    for i = 1:binnum
        string = ['numbins = ' num2str(binnum) '; binnum #: ' num2str(i)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        beginboxcar = (((i*interval)-interval)+1);
        endboxcar = (((i*interval)-interval)+binlength);
        induse = ind(beginboxcar:endboxcar);
        box_corrmat = corrcoef(tc_concat{s}(induse,:));
        diff_corrmat = box_corrmat - good_corrmat;
        mean_diff(i) = mean(diff_corrmat(mask));
        mean_FD_box(i) = mean(FD_concat{s}(induse));%mean(DV_concat{s}(induse));
        
        for p = 1:250
            randind = randperm(length(gooddata_induse));%randperm(size(tc_concat{s},1));
            tc_concat_rand = tc_concat{s}(gooddata_induse(randind(1:binlength)),:);
            rand_corrmat = corrcoef(tc_concat_rand);
            diff_rand_corrmat = rand_corrmat - good_corrmat;
            mean_rand_diff(p,i) = mean(diff_rand_corrmat(mask));
        end
        
        
    end
    disp(' ')
    
    
    
    %% Display effect of FD on delta r
    h = figure('Position',[1994 324 1155 746],'Color','white');
    mean_FD_box_rand = repmat(mean_FD_box,[250 1]);
    
    plot(mean_FD_box_rand,mean_rand_diff,'.k','MarkerSize',10)
    
    hold
    
    for i = 1:binnum
        
        cross_subject_mean_FD_box = [cross_subject_mean_FD_box; mean_FD_box(i)];
        cross_subject_mean_diff_v_rand = [cross_subject_mean_diff_v_rand; (nnz(mean_diff(i)>mean_rand_diff(:,i)) / size(mean_rand_diff,1))];
        
        if nnz(mean_diff(i)<mean_rand_diff(:,i))<13
            plot(mean_FD_box(i),mean_diff(i),'.r','MarkerSize',20)
        else
            plot(mean_FD_box(i),mean_diff(i),'.g','MarkerSize',10)
        end
    end
    ylim([-.1 .1])
    ylabel('delta','FontWeight','bold','Fontsize',14)
    xlabel('FD','FontWeight','bold','Fontsize',14)
    set(gca,'FontWeight','bold','Fontsize',14)
    
    drawnow
    
    
%     %% Look for distance dependent artifact
%     maskmat = ones(length(values));
%     maskmat = triu(maskmat,1);
%     
%     lower_FDthresh = [0 FDthresh(1:end-1)];
%     
%     for f = 1:length(FDthresh)
%         
%         induse = (FD_concat{s}<FDthresh(f)); %& (FD_concat{s}>lower_FDthresh(f));
%         corrmat_scrub = corrcoef(tc_concat{s}(induse,:));
%         
%         %Calculate delta r
%         diffmat = all_corrmat-corrmat_scrub;
%         
%         %F0 = fit(coords_dist(logical(maskmat)),diffmat(logical(maskmat)),'poly1');
%         %B = regress(diffmat(logical(maskmat)),[coords_dist(logical(maskmat)) ones(nnz(maskmat),1)]);
%         B = regress(diffmat(coords_dist_short & logical(maskmat)),[coords_dist(coords_dist_short & logical(maskmat)) ones(nnz(coords_dist_short & maskmat),1)]);
%         allslopes(s,f) = B(1);
        
        
%         yhat = B(1) .* coords_dist(logical(maskmat)) + B(2);
%         
%         %Sort r
%         distmat_vec = coords_dist(logical(maskmat));
%         [distmat_vec_sort ind] = sort(distmat_vec,'ascend');
%         diffmat_vec = diffmat(logical(maskmat));
%         diffmat_vec_sort = diffmat_vec(ind);
%         
%         finedistbins=[10:2:170];
%         finedistbincenters=[11:2:169];
%         
%         plotarray = calcdistrmat(diffmat,finedistbins,coords_dist);
%         
%         clear distbin_center diffbin_val
%         bins = 1:300:34716;
%         for b = 1:length(bins)-1
%             distbin_center(b) = mean(distmat_vec_sort(bins(b):bins(b+1)));
%             diffbin_val(b) = mean(diffmat_vec_sort(bins(b):bins(b+1)));
%         end
%         
%         
%         h = figure('Color','white','Position',[2147 256 916 660]);
%         hold
%         
%         plot(coords_dist(logical(maskmat)),diffmat(logical(maskmat)),'.')
%         %plot(distbin_center,diffbin_val,'-r','LineWidth',2)
%         plot(finedistbincenters,plotarray,'-r','LineWidth',3)
%         %hline_new(0,'k',2)
%         plot(coords_dist(logical(maskmat)),yhat,'-g')
%         xlabel('distance (mm)','FontWeight','bold','FontSize',14)
%         ylabel('delta r','FontWeight','bold','FontSize',14)
%         title(['FD thresh = ' num2str(FDthresh(f))],'FontWeight','bold','FontSize',14)
%         set(gca,'FontWeight','bold','FontSize',14)
%         
%         waitforbuttonpress
%         close(h)
%    end
    %figure;plot(FDthresh,slopes,'k')
    
    
    
end

%figure; plot(FDthresh,allslopes([1:5 7:11],:))
        
    
%figure;plot(cross_subject_mean_FD_box,cross_subject_mean_diff_v_rand,'.')
    
    
    





























