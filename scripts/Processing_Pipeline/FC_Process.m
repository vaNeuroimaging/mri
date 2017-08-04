function FC_Process(subject,BOLDlist,FDthresh,TR,outfolder,outprefix)
%FC_Process(subject,BOLDlist,FDthresh,TR,outfolder,outprefix)

brainradius = 50;
contiguoustimepoints = 5;
lowpassfilt = .08; %Hz
highpassfilt = .009; %Hz
minpctframes = .1;





%Calculate filter properties
lopasscutoff=lowpassfilt/(0.5/TR); 
hipasscutoff=highpassfilt/(0.5/TR);
[butta, buttb]=butter(1,[hipasscutoff lopasscutoff]);



voxels_toprocess = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_brainmask_dil1_mask_333.nii.gz']); voxels_toprocess = logical(voxels_toprocess.img);
%voxels_toprocess = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_brainmask_mask_333.nii.gz']); voxels_toprocess = logical(voxels_toprocess.img);

%Get nuisance regressor masks
wholebrainmask = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_brainmask_mask_333.nii.gz']); wholebrainmask = logical(wholebrainmask.img(voxels_toprocess));

WMero=4;
done = 0;
while done==0
    WMmask = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_cerebralwm_ero' num2str(WMero) '_mask_333.nii.gz']); WMmask = logical(WMmask.img(voxels_toprocess));
    if nnz(WMmask)==0;
        WMero = WMero-1;
    else
        done = 1;
    end
end

CSFero=1;
done = 0;
while done==0
    CSFmask = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_cerebralwm_ero' num2str(CSFero) '_mask_333.nii.gz']); CSFmask = logical(CSFmask.img(voxels_toprocess));
    if nnz(CSFmask)==0;
        CSFero = CSFero-1;
    else
        done = 1;
    end
end




%Get list of files and motion params
[runs,sessions,motionparamsfiles] = textread(BOLDlist,'%s%s%s');
nruns = length(motionparamsfiles);

all_tmask = [];
FD = [];
run_tmasks = cell(length(motionparamsfiles),0);
runIDs = [];
sessionIDs = [];
motion_regressors = [];

for runnum = 1:nruns
    
    
    %Calculate motuion regressors
    thisrun_params = load(motionparamsfiles{runnum});
    
    thisrun_rot = thisrun_params(:,1:3);
    thisrun_rot_mm = thisrun_rot * brainradius;
    thisrun_params(:,1:3) = thisrun_rot_mm;
    thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];
    
    %demean and detrend motion regressors
    thisrun_params_demeandetrend = demean_detrend(thisrun_params');
    thisrun_params_delta_demeandetrend = demean_detrend(thisrun_params_delta');
    motion_regressors = [motion_regressors [thisrun_params_demeandetrend ; thisrun_params_delta_demeandetrend; thisrun_params_demeandetrend.^2]];
    
    
    %Calculate FD and tmask
    thisrun_FD{runnum} = [sum(abs(thisrun_params_delta),2)]';
    FD = [FD ; thisrun_FD{runnum}'];
    thisrun_tmask = (thisrun_FD{runnum}') < FDthresh;
    
    %Apply contiguous frame criterion
    chunkedges = [1; diff([thisrun_tmask;0])];
    chunkstarts = find(chunkedges==1); chunkstops = (find(chunkedges==-1) -1);
    chunksizes = zeros(length(chunkstarts),1);
    chunkID = zeros(size(thisrun_tmask));
    for chunknum = 1:length(chunkstarts)
        chunkID(chunkstarts(chunknum) : chunkstops(chunknum)) = chunknum;
        chunksizes(chunknum) = nnz(chunkID==chunknum);
        if chunksizes(chunknum) < contiguoustimepoints
            thisrun_tmask(chunkID==chunknum) = 0;
        end
    end
    
    if (nnz(thisrun_tmask) / numel(thisrun_tmask)) < minpctframes;
        thisrun_tmask(:) = 0;
    end
    
    run_tmasks{runnum} = thisrun_tmask;
    all_tmask = [all_tmask ; thisrun_tmask];
    
    
    %Load data
    thisrun_data = load_untouch_nii_2D(runs{runnum});
    if runnum==1
        output_structure = thisrun_data; output_structure.img = [];
    end
    thisrun_data = thisrun_data.img(voxels_toprocess,:);
    
    %Demean and detrend
    if any(run_tmasks{runnum})
        thisrun_data = demean_detrend(thisrun_data,run_tmasks{runnum});
    end
    
    %Save in one big variable
    if runnum==1
        data = thisrun_data;
    else
        data = [data thisrun_data];
    end
    
    %Track run IDs
    runIDs = [runIDs; repmat(runnum,size(thisrun_data,2),1)];
    sessionIDs = [sessionIDs; repmat(num2str(sessions{runnum}),size(thisrun_data,2),1)];
    
    clear thisrun_data
    
end


dlmwrite([outfolder '/' outprefix '_all_tmask.txt'],all_tmask,'delimiter',' ')
dlmwrite([outfolder '/' outprefix '_all_FD.txt'],FD,'delimiter',' ')
dlmwrite([outfolder '/' outprefix '_runs_sessions.txt'],[runIDs str2num(sessionIDs)],'delimiter',' ')



%Get physiological regressors
WBtimecourse = mean(data(wholebrainmask,:),1);
WMtimecourse = mean(data(WMmask,:),1);
CSFtimecourse = mean(data(CSFmask,:),1);
physio_regressors = [WBtimecourse ; WMtimecourse ; CSFtimecourse];
physio_regressors_delta = diff([zeros(3,1) physio_regressors],1,2);

%Add to motion regressors
all_regressors = [motion_regressors ; physio_regressors ; physio_regressors_delta]';



%Regress out nuisance signals
data = regress_nuisance(data,all_regressors,all_tmask);



for runnum = 1:nruns
    
    thisrun_data = data(:,runIDs==runnum);
        
    if any(run_tmasks{runnum})
        
        %Interpolate data across censored frames in preparation for filtering
        ofac=8;
        hifac=1;
        TRtimes=([1:size(thisrun_data,2)]').*TR;
        
        if numel(TRtimes)<150
            voxbinsize=5000;
        elseif (numel(TRtimes)>=150 && numel(TRtimes)<500)
            voxbinsize=500;
        elseif numel(TRtimes)>=500
            voxbinsize=50;
        end
        voxbin=[1:voxbinsize:size(thisrun_data,1) size(thisrun_data,1)];
        
        interpolated=zeros(size(thisrun_data,2),size(thisrun_data,1));
        
        % gotta bin by voxels: 5K is ~15GB, 15K is ~40GB at standard
        % run lengths. 5K is ~15% slower but saves 2/3 RAM, so that's
        % the call.
        for v=1:numel(voxbin)-1 % this takes huge RAM if all voxels
            interpolated(:,voxbin(v):voxbin(v+1))=Interpolate_intmask(TRtimes(run_tmasks{runnum}),thisrun_data(voxbin(v):voxbin(v+1),run_tmasks{runnum})',TRtimes,TR,ofac,hifac);
        end
        
        interpolated=interpolated';
        
        thisrun_data(:,~run_tmasks{runnum})=interpolated(:,~run_tmasks{runnum});
        clear interpolated
    end
    
    %Filter data
    nvox = size(thisrun_data,1);
    padsize = 1000;
    
    padded_data = [zeros(padsize,nvox) ; thisrun_data' ; zeros(padsize,nvox)];
    filtered = filtfilt(butta,buttb,double(padded_data));
    thisrun_data = filtered(padsize+1 : end-padsize,:)';
    
    
    %Demean and detrend data once more
    thisrun_data = demean_detrend(thisrun_data,run_tmasks{runnum});
    
    data(:,runIDs==runnum) = thisrun_data;
    
    
    
end

output_structure.hdr.dime.dim(5) = size(data,2);
output_structure.img = zeros(prod(output_structure.hdr.dime.dim(2:4)),size(data,2));
output_structure.img(voxels_toprocess,:) = data; 
clear data
outfile = [outfolder '/' outprefix '_fc_processed.nii.gz'];
save_untouch_nii_2D(output_structure,outfile);

output_structure.img = output_structure.img(:,logical(all_tmask));
output_structure.hdr.dime.dim(5) = nnz(all_tmask);
outfile = [outfolder '/' outprefix '_fc_processed_tmasked.nii.gz'];
save_untouch_nii_2D(output_structure,outfile);


end





function [bold,tempbetas] = demean_detrend(bold,tmask)
if ~exist('tmask')
    tmask = true(1,size(bold,2));
end
[vox,ts] = size(bold);
linreg=[repmat(1,[ts 1]) linspace(0,1,ts)'];
tempboldcell=num2cell(bold(:,logical(tmask))',1);
linregcell=repmat({linreg(logical(tmask),:)},[1 vox]);
tempbetas = cellfun(@mldivide,linregcell,tempboldcell,'uniformoutput',0);
tempbetas=cell2mat(tempbetas);
tempbetas=tempbetas';
tempintvals=tempbetas*linreg';
bold=bold-tempintvals;

end



function [tempimg, zb, newregs] = regress_nuisance(tempimg,totregs,tot_tmask)

[vox ts]=size(tempimg);
zlinreg=totregs(logical(tot_tmask),:); % only desired data
[zlinreg DMDTB]=demean_detrend(zlinreg'); % obtain fits for desired data
zlinreg=zlinreg';
zstd=std(zlinreg); % calculate std
zmean=mean(zlinreg);
zlinreg=zlinreg-(repmat(zmean,[size(zlinreg,1) 1]))./(repmat(zstd,[size(zlinreg,1) 1])); % zscore

linreg=[repmat(1,[ts 1]) linspace(0,1,ts)'];
newregs=DMDTB*linreg'; % predicted all regressors demean/detrend
newregs=totregs-newregs'; % these are the demeaned detrended regressors
newregs=newregs-(repmat(zmean,[size(newregs,1) 1]))./(repmat(zstd,[size(newregs,1) 1])); % zscore

% now we have z-scored, detrended good and all regressors.

% demean and detrend the desired data
zmdtimg=tempimg(:,logical(tot_tmask));
[zmdtimg zmdtbetas]=demean_detrend(zmdtimg);

% calculate betas on the good data
tempboldcell=num2cell(zmdtimg',1);
zlinregcell=repmat({zlinreg},[1 vox]);
zb = cellfun(@mldivide,zlinregcell,tempboldcell,'uniformoutput',0);
zb=cell2mat(zb);

% demean and detrend all data using good fits
[zmdttotimg]=zmdtbetas*linreg';
zmdttotimg=tempimg-zmdttotimg;

% calculate residuals on all the data
zb=zb';
tempintvals=zb*newregs';
tempimg=zmdttotimg-tempintvals;

end


function [H,f,s,c,tau,w] = Interpolate_intmask(t,h,TH,Tr,ofac,hifac)

%Input t is a column vector listing the time points for which observations
%are present.  Input h is a matrix with observations in columns and the
%number of rows equals the number the time points.  For our purposes number
%of voxels = number of columns.  Ofac = oversampling frequency (generally
%>=4), hifac = highest frequency allowed.  hifac = 1 means 1*nyquist limit
%is highest frequency sampled.  
%Lasted edited:  Anish Mitra, October 25 2012

N = size(h,1); %Number of time points
T = max(t) - min(t); %Total time span

%calculate sampling frequencies
f = (1/(T*ofac):1/(T*ofac):hifac*N/(2*T)).';

%angular frequencies and constant offsets
w = 2*pi*f;
tau = atan2(sum(sin(2*w*t.'),2),sum(cos(2*w*t.'),2))./(2*w);

%spectral power sin and cosine terms
cterm = cos(w*t.' - repmat(w.*tau,1,length(t)));
sterm = sin(w*t.' - repmat(w.*tau,1,length(t)));

num_tseries = size(h,2); %Number of time series

% D = bsxfun(@minus,h,mean(h));  %This line involved in normalization
D = h;
D = reshape(D,1,N,num_tseries);


%% C_final = (sum(Cmult,2).^2)./sum(Cterm.^2,2);
% This calculation is done by separately for the numerator, denominator,
% and the division

Cmult = bsxfun(@times, cterm,D);
%rewrite the above line with bsxfun to optimize further?
%bsxfun(@power,sum(Cmult,2),2) = sum(Cmult.^2,2) = numerator
% numerator = bsxfun(@power,sum(Cmult,2),2);
numerator = sum(Cmult,2); %Modify the numerator to get the cw term

%sum(bsxfun(@power,Cterm,2),2) = sum(Cterm.^2,2) = denominator
denominator = sum(bsxfun(@power,cterm,2),2); %use Cterm in place of cterm to make it exactly the denominator in the original expression
C_final_new = bsxfun(@rdivide,numerator,denominator);
c = C_final_new;
clear numerator denominator cterm Cmult C_final_new

%% Repeat the above for Sine term
Smult = bsxfun(@times, sterm,D);
% S_final = (sum(Smult,2).^2)./sum(Sterm.^2,2);
% numerator = bsxfun(@power,sum(Smult,2),2);
numerator = sum(Smult,2); %Modify the numerator to get the sw term
denominator = sum(bsxfun(@power,sterm,2),2);
S_final_new = bsxfun(@rdivide,numerator,denominator);
s = S_final_new;
clear numerator denominator sterm Smult S_final_new

% %% Power = C_final + S_final; 
% Power = C_final_new + S_final_new;
% Power_reshaped = reshape(Power,size(Power,1),num_tseries);
% Power_final = bsxfun(@rdivide,Power_reshaped,2*var(h)); %Normalize the power
% clearvars -except Power_final f


% The inverse function to re-construct the original time series
Time = TH';
T_rep = repmat(Time,[size(f,1),1,size(h,2)]);
% T_rep = bsxfun(@minus,T_rep,tau);
% s(200) = s(199);
w = 2*pi*f;
prod = bsxfun(@times,T_rep,w);
sin_t = sin(prod);
cos_t = cos(prod);
sw_p = bsxfun(@times,sin_t,s);
cw_p = bsxfun(@times,cos_t,c);
S = sum(sw_p);
C = sum(cw_p);
H = C + S;
H = reshape(H,size(Time,2),size(h,2));

%Normalize the reconstructed spectrum, needed when ofac > 1
Std_H = std(H);
Std_h = std(h);
norm_fac = Std_H./Std_h;
H = bsxfun(@rdivide,H,norm_fac);

end
