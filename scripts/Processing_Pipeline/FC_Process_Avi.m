function FC_Process_Avi(subject,BOLDlist,FDthresh,TR,outfolder,outprefix)
%FC_Process(subject,BOLDlist,FDthresh,TR,outfolder,outprefix)
%
%Inputs:
% subject - text string specifying name of subject. Will be used to find
%  the appropriate nuisance masks.
% BOLDlist - a text file with the following format:
%       Run1_preprocesseddata.nii.gz  Run1_sessionnumber  Run1_motionparams.par
%       Run2_preprocesseddata.nii.gz  Run2_sessionnumber  Run2_motionparams.par
%       Run3_preprocesseddata.nii.gz  Run3_sessionnumber  Run3_motionparams.par
%       ...etc
% FDthresh - the FD cutoff to use for censoring, in mm
% TR - the TR of the data, in seconds.
% outfolder - location where the data will be written
% outprefix - string to append to beginning of output file

brainradius = 50;
contiguoustimepoints = 5; %to pass the censoring threshold, [contiguoustimepoints] consecutive frames must be below the threshold. 
lowpassfilt = .08; %Hz
highpassfilt = .009; %Hz
minpctframes = .1; %Runs will not be included if less than [minpctframes] of the frames survive thresholding 
CSF_PCA_cutoff = .75; %Percent variance in CSF voxel timecourses explained by to-be-included CSF PCA components 


voxels_toprocess = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_brainmask_dil1_mask_333.nii.gz']); voxels_toprocess = logical(voxels_toprocess.img);

%Get nuisance regressor masks
wholebrainmask = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_brainmask_mask_333.nii.gz']); wholebrainmask = logical(wholebrainmask.img(voxels_toprocess));
WMmask = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_cerebralwm_mask_333_6mm_fromGMCSF.nii.gz']); WMmask = logical(WMmask.img(voxels_toprocess));
CSFmask = load_untouch_nii_2D(['/home/data/subjects/' subject '/freesurfer/nusmask/aparc+aseg_CSF_mask_333_6mm_fromGMWM.nii.gz']); CSFmask = logical(CSFmask.img(voxels_toprocess));



%------------------------------------------------



%Calculate filter properties
lopasscutoff=lowpassfilt/(0.5/TR); 
hipasscutoff=highpassfilt/(0.5/TR);
[butta, buttb]=butter(1,[hipasscutoff lopasscutoff]);



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
%     thisrun_params = load(motionparamsfiles{runnum});
%     thisrun_rot = thisrun_params(:,1:3);
%     thisrun_rot_mm = thisrun_rot * brainradius;
%     thisrun_params(:,1:3) = thisrun_rot_mm;
    
    thisrun_params = load(motionparamsfiles{runnum});
    thisrun_rot = thisrun_params(:,1:3);
    thisrun_rot_mm = thisrun_rot * brainradius;
    thisrun_params(:,1:3) = thisrun_rot_mm;
    %figure;
    %plot(thisrun_params_unfilt)
    
    thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];
    
    %demean and detrend motion regressors
    thisrun_params_demeandetrend = demean_detrend(thisrun_params');
    thisrun_params_delta_demeandetrend = demean_detrend(thisrun_params_delta');
    thisrun_motion_regressors = [thisrun_params_demeandetrend ; thisrun_params_delta_demeandetrend; thisrun_params_demeandetrend.^2];
    
    
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
    
    %If the run has data
    if any(run_tmasks{runnum})
        %Demean and detrend
        thisrun_data = demean_detrend(thisrun_data,run_tmasks{runnum});
        thisrun_motion_regressors = demean_detrend(thisrun_motion_regressors,run_tmasks{runnum});
        
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
        
        
        interpolated=Interpolate_intmask(TRtimes(run_tmasks{runnum}),thisrun_motion_regressors(:,run_tmasks{runnum})',TRtimes,TR,ofac,hifac);
        interpolated=interpolated';
        thisrun_motion_regressors(:,~run_tmasks{runnum})=interpolated(:,~run_tmasks{runnum});
        clear interpolated
        
        
    else
        thisrun_data = demean_detrend(thisrun_data,ones(size(run_tmasks{runnum})));
        thisrun_motion_regressors = demean_detrend(thisrun_motion_regressors,ones(size(run_tmasks{runnum})));
    end
    
    %Filter data
    nvox = size(thisrun_data,1);
    padsize = 1000;
    padded_data = [zeros(padsize,nvox) ; thisrun_data' ; zeros(padsize,nvox)];
    filtered = filtfilt(butta,buttb,double(padded_data));
    thisrun_data = filtered(padsize+1 : end-padsize,:)';
    
    %Filter regressors
    nregress = size(thisrun_motion_regressors,1);
    padsize = 1000;
    padded_data = [zeros(padsize,nregress) ; thisrun_motion_regressors' ; zeros(padsize,nregress)];
    filtered = filtfilt(butta,buttb,double(padded_data));
    thisrun_motion_regressors = filtered(padsize+1 : end-padsize,:)';
    
    
    %Demean and detrend data once more
    if any(run_tmasks{runnum})
        thisrun_data = demean_detrend(thisrun_data,run_tmasks{runnum});
        thisrun_motion_regressors = demean_detrend(thisrun_motion_regressors,run_tmasks{runnum});
    end
    
    
    
    %Save in one big variable
    if runnum==1
        data = thisrun_data;
        motion_regressors = thisrun_motion_regressors;
    else
        data = [data thisrun_data];
        motion_regressors = [motion_regressors thisrun_motion_regressors];
    end
    
    
    %Track run IDs
    runIDs = [runIDs; repmat(runnum,size(thisrun_data,2),1)];
    sessionIDs = [sessionIDs; repmat(str2num(sessions{runnum}),size(thisrun_data,2),1)];
    
    clear thisrun_data
    
end

data(isnan(data)) = 0;

dlmwrite([outfolder '/' outprefix '_all_tmask.txt'],all_tmask,'delimiter',' ')
dlmwrite([outfolder '/' outprefix '_all_FD.txt'],FD,'delimiter',' ')
dlmwrite([outfolder '/' outprefix '_runs_sessions.txt'],[runIDs sessionIDs],'delimiter',' ')


%Get physiological regressors
WBtimecourse = mean(data(wholebrainmask,:),1);
WMtimecourse = mean(data(WMmask,:),1);
[PCs,PCscore,eigvals_sort,~,pct_explained] = pca(data(CSFmask,:)');
cum_eigvals_perc = cumsum(pct_explained);
numcomps_toinclude = max(find(cum_eigvals_perc < (CSF_PCA_cutoff*100)))+1;
CSFregressors = PCscore(:,1:numcomps_toinclude)';
physio_regressors = [WBtimecourse ; WMtimecourse ; CSFregressors];
physio_regressors_delta = [zeros(size(physio_regressors,1),1) diff(physio_regressors,1,2)];
physio_regressors_plusdeltas = [physio_regressors; physio_regressors_delta];

%Add to motion regressors
all_regressors = [motion_regressors ; physio_regressors_plusdeltas]';


%Regress out nuisance signals
data = regress_nuisance(data,all_regressors,all_tmask);


%Save outputs
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


string = [subject ', ' outprefix ': ' num2str(nnz(all_tmask) / numel(all_tmask)) '% (' num2str(nnz(all_tmask)) ' of ' num2str(numel(all_tmask)) ') frames retained'];
%disp(string)
figure;
plot(FD,'k')
hold on
plot((ones(length(FD),1) .* FDthresh),'--r')
nantmask = single(~all_tmask); nantmask(logical(nantmask)) = NaN;
plot(nantmask+.001,'b','Linewidth',5)
title(string)
export_fig(gca,[outfolder '/' outprefix '_all_FD.pdf'])


end





function [bold,tempbetas] = demean_detrend(bold,tmask)
if ~exist('tmask')
    tmask = true(1,size(bold,2));
end
[vox,ts] = size(bold);
linreg=[ones([ts 1]) linspace(0,1,ts)'];
tempboldcell=num2cell(bold(:,logical(tmask))',1);
linregcell=repmat({linreg(logical(tmask),:)},[1 vox]);
tempbetas = cellfun(@mldivide,linregcell,tempboldcell,'uniformoutput',0);
tempbetas=cell2mat(tempbetas);
tempbetas=tempbetas';
tempintvals=tempbetas*linreg';
bold=bold-tempintvals;

end



function [tempimg, zb, newregs] = regress_nuisance(tempimg,totregs,tot_tmask)

[vox, ts]=size(tempimg);
allzeros = logical(all(totregs(logical(tot_tmask),:)==0,1));
totregs(:,allzeros) = [];
zlinreg=totregs(logical(tot_tmask),:); % only desired data
[zlinreg, DMDTB]=demean_detrend(zlinreg'); % obtain fits for desired data
zlinreg=zlinreg';
zstd=std(zlinreg); % calculate std
zmean=mean(zlinreg);
zlinreg=zlinreg-(repmat(zmean,[size(zlinreg,1) 1]))./(repmat(zstd,[size(zlinreg,1) 1])); % zscore

linreg=[ones([ts 1]) linspace(0,1,ts)'];
newregs=DMDTB*linreg'; % predicted all regressors demean/detrend
newregs=totregs-newregs'; % these are the demeaned detrended regressors
newregs=newregs-(repmat(zmean,[size(newregs,1) 1]))./(repmat(zstd,[size(newregs,1) 1])); % zscore

% now we have z-scored, detrended good and all regressors.

% demean and detrend the desired data
zmdtimg=tempimg(:,logical(tot_tmask));
[zmdtimg, zmdtbetas]=demean_detrend(zmdtimg);

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
