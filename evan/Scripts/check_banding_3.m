subjects = textread('/home/data/subjects/processing_list_10min.txt','%s');

for s = length(subjects)
%subjects{1} = 'MAV033';
%s=1;
    prevstring = [];
    %cd(['/home/data/subjects/' subjects{s} '/fc_processed/'])
    cd(['/home/data/subjects/' subjects{s} '/preprocessed/'])
    
    %data = load_untouch_nii('restingstate_fc_processed_tmasked.nii.gz');
    data = load_untouch_nii('restingstate_1_mcf.nii.gz');
    
    data = data.img;
    
    %wm = load_untouch_nii('../freesurfer/nusmask/aparc+aseg_cerebralwm_ero0_mask_333.nii.gz');
    wm = load_untouch_nii('aparc+aseg_cerebralwm_ero0_mask_restingstate1space.nii.gz');
    wm = wm.img;
    wm(wm==0) = NaN;
    data = data .* repmat(wm,[1 1 1 size(data,4)]);
    
    mean_power{s} = zeros(size(data,4),1);
    mean_power_pct{s} = zeros(size(data,4),1);
    for t = 1:size(data,4);
        string = [subjects{s} ', timepoint ' num2str(t)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        timepointdata = data(:,:,:,t);
        zvec = nanmean(timepointdata,2);
        zvec = squeeze(nanmean(zvec,1));
        zvec(isnan(zvec)) = [];
        

            
            
            

        
        mean_power{s}(t) = bandpower(zvec,1,[1/10 1/3]);
        mean_power_pct{s}(t) = mean_power{s}(t) / bandpower(zvec);
        
       %         figure
        %         plotSpread(bin_rsq,'distributionColors',{'k'},'MarkerSize',15)
        %     title(['Timepoint ' num2str(t) ': median r-sq=' num2str(median_rsqs(t))])
        %
        %     waitforbuttonpress
        %
        %     close all
    end
    disp(' ')
    
    %out(:,:,1:55,:) = tempout;
end

