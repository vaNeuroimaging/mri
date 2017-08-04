%subjects = textread('/home/data/subjects/processing_list_10min.txt','%s');

%for s = 3:length(subjects)
subjects{1} = 'MAV033';
s=1;
    prevstring = [];
    cd(['/home/data/subjects/' subjects{s} '/fc_processed/'])
    
    data = load_untouch_nii('restingstate_fc_processed_tmasked.nii.gz');
    
    data = data.img;
    
    wm = load_untouch_nii('../freesurfer/nusmask/aparc+aseg_cerebralwm_ero0_mask_333.nii.gz');
    wm = wm.img;
    
    ybins = [20:1:40];
    %ybinsize = 5;
    
    out = zeros(size(data));
    
    data = data(:,:,1:55,:) .* repmat(wm(:,:,1:55,:),[1 1 1 size(data,4)]);
    tempout = zeros(size(data));
    median_rsqs{s} = zeros(size(data,4),1);
    for t = 1:30%size(data,4);
        string = [subjects{s} ', timepoint ' num2str(t)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        for binnum = 1:(length(ybins)-1)
            %timepointdata = data(:,:,:,t) .* wm;
            %timepointdata = timepointdata(:,30:40,30:end);
            timepointdata = data(:,ybins(binnum):(ybins(binnum+1)-1),:,t);
            datainds = find(timepointdata~=0);
            if length(datainds) > 10
            [~,~,z] = ind2sub(size(timepointdata),datainds);
            inbraindata = single(timepointdata(datainds));
            
            [FO, G] = fit(z,inbraindata,'sin1');
            resid = inbraindata - FO(z);
            
            timepointdata(datainds) = resid;
            tempout(:,ybins(binnum):(ybins(binnum+1)-1),:,t) = timepointdata;
            
            bin_rsq(binnum,1) = G.rsquare;
            end
        end
        
        %         timepointdata = data(:,:,:,t) .* wm;
        %         timepointdata = timepointdata(:,30:40,30:end);
        %         datainds = find(timepointdata~=0);
        %         [~,~,z] = ind2sub(size(timepointdata),datainds);
        %         inbraindata = timepointdata(datainds);
        %
        %         [FO, G] = fit(z,inbraindata,'sin1');
        %
        %         uniquez = unique(z);
        %         fitted = FO(uniquez);
        %
        %         plot(z,inbraindata,'k.')
        %         hold on
        %         plot(uniquez,fitted,'b-')
        %         title(['Timepoint ' num2str(t) ': r-sq=' num2str(G.rsquare)])
        median_rsqs{s}(t) = median(bin_rsq);
        
        %         figure
        %         plotSpread(bin_rsq,'distributionColors',{'k'},'MarkerSize',15)
        %     title(['Timepoint ' num2str(t) ': median r-sq=' num2str(median_rsqs(t))])
        %
        %     waitforbuttonpress
        %
        %     close all
    end
    disp(' ')
    
    out(:,:,1:55,:) = tempout;
%end

