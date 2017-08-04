function stripecheck(directory)

warning off


FDthresh = .2;
brainradius = 50;
contiguoustimepoints = 5;


FD = [];
newFD = [];

localdir = [directory '/temp/'];
mkdir(localdir)
system(['cp ' directory '/restingstate*.nii.gz ' localdir])

files = dir([localdir '/*.nii.gz']);
%files = struct2cell(files); files = files(1,:);

%copyfile([localdir '/*.nii.gz'],localdir)
%ign = rmdir(tempdir,'s');

%thisseq_files = dir([localdir '/*.nii.gz']);
%disp(['Converted ' num2str(length(files)) ' new runs and found ' num2str(length(thisseq_files)-length(files)) ' previous runs in ' localdir '. Getting motion estimates.'])

string = 'fslview ';

for f = 1:length(files);
    
    data = [localdir files(f).name(1:end-7)];
    
    paramsfile = [data '_mcf.par'];
    
    if ~exist(paramsfile,'file')
        [fail, result] = system(['mcflirt -in ' data ' -refvol 0 -plots; '...
            'flirt -in ' data '_mcf -ref /usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz -dof 12 -omat BOLD2MNI.mat; ' ...
            'flirt -in ' data '_mcf -ref /usr/local/fsl/data/standard/MNI152lin_T1_1mm_brain.nii.gz -applyisoxfm 3 -init BOLD2MNI.mat -out ' data '_mcf_MNI;']);
        if fail
            disp(['WARNING: could not process ' data ': ' result])
        end
    end
    
    a = load_untouch_nii_2D([data '_mcf_MNI.nii.gz']);
    b = a;
    b.img = demean_detrend(a.img);
    save_untouch_nii_2D(b,[data '_mcf_MNI_dmdt.nii.gz'])
    string = [string data '_mcf_MNI_dmdt -l render3 -b -20,20 '];
    %system(['fslview ' data '_mcf_MNI_dmdt -l render3 -b -20,20 & '])
    
    
%     FC_Process_noseg([data '_mcf_MNI.nii.gz'],[data '_mcf.par']);
%     
%     
%     
%     data = load_untouch_nii([data '_mcf_MNI_fc_processed_tmasked.nii.gz']);
%     
%     data = data.img;
%     
%     wm = load_untouch_nii('/Users/user/Documents/MATLAB/avg152T1_white_thr240_333.nii.gz');
%     wm = wm.img;
%     wm(wm==0) = NaN;
%     data = data .* repmat(wm,[1 1 1 size(data,4)]);
%     
%     mean_power = zeros(size(data,4),1);
%     mean_power_pct = zeros(size(data,4),1);
%     for t = 1:size(data,4);
%         string = [subjects{s} ', timepoint ' num2str(t)];
%         fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%         prevstring = string;
%         
%         timepointdata = data(:,:,:,t);
%         zvec = nanmean(timepointdata,2);
%         zvec = squeeze(nanmean(zvec,1));
%         zvec(isnan(zvec)) = [];
%         
%         mean_power(t) = bandpower(zvec,1,[1/10 1/3]);
%         mean_power_pct(t) = mean_power{s}(t) / bandpower(zvec);
%         
%     end
%     disp(' ')
%     
%     figure;
%     plot(mean_power)
    
    
    
    
end

system(string)







