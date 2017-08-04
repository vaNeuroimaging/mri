subjects = textread('/home/data/subjects/processing_list_10min.txt','%s');

for s = 1:length(subjects)
%subjects{1} = 'MAV033';
%s=1;
    prevstring = [];
    cd(['/home/data/subjects/' subjects{s} '/fc_processed/'])
    
    data = load_untouch_nii('restingstate_fc_processed_tmasked.nii.gz');
    
    data = data.img;
    
    wm = load_untouch_nii('../freesurfer/nusmask/aparc+aseg_cerebralwm_ero0_mask_333.nii.gz');
    wm = wm.img;
    wm(wm==0) = NaN;
    data = data .* repmat(wm,[1 1 1 size(data,4)]);
    
    ybins = [10:2:55];
    
    mean_power{s} = zeros(size(data,4),1);
    mean_power_pct{s} = zeros(size(data,4),1);
    for t = 1:size(data,4);
        string = [subjects{s} ', timepoint ' num2str(t)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        for binnum = 1:(length(ybins)-1)
            
            timepointdata = data(:,ybins(binnum):(ybins(binnum+1)-1),:,t);
            
            zvec = nanmean(timepointdata,2);
            zvec = squeeze(nanmean(zvec,1));
            zvec(isnan(zvec)) = [];
            
            powerestimate(binnum) = bandpower(zvec,1,[1/10 1/3]);
            pct_powerestimate(binnum) = powerestimate(binnum) / bandpower(zvec);
            
        end
        
        mean_power{s}(t) = mean(powerestimate);
        mean_power_pct{s}(t) = mean(pct_powerestimate);
    end
    disp(' ')
end

