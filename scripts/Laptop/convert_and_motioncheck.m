function convert_and_motioncheck(directory)

warning off

fsldir = '/usr/local/fsl/'; 
dcm2niidir = '/Users/user/Documents/mricron/';
localdir = '/Users/user/Documents/motioncheck/';

system('export FSLDIR=/usr/local/fsl/; export PATH=/usr/local/fsl/bin/:$PATH; source /usr/local/fsl/etc/fslconf/fsl.sh')

addpath(fsldir);
addpath(dcm2niidir);

disp('Converting...')
tempdir = [localdir 'temp/'];
mkdir(tempdir);
[failed, message] = system([dcm2niidir 'dcm2nii -o ' tempdir ' ' directory]);
if failed
    disp('Dicom conversion FAILED!')
    disp(message)
end



FDthresh = .2;
brainradius = 50;
contiguoustimepoints = 5;


FD = [];
newFD = [];

newfiles = dir([tempdir '/*.nii.gz']);
newfiles = struct2cell(newfiles); newfiles = newfiles(1,:);

copyfile([tempdir '/*.nii.gz'],localdir)
ign = rmdir(tempdir,'s');

thisseq_files = dir([localdir '/*.nii.gz']);
disp(['Converted ' num2str(length(newfiles)) ' new runs and found ' num2str(length(thisseq_files)-length(newfiles)) ' previous runs in ' localdir '. Getting motion estimates.'])

for f = 1:length(thisseq_files);
    
    data = [localdir thisseq_files(f).name];
    
    paramsfile = [data(1:end-7) '_mcf.par'];
    
    if ~exist(paramsfile,'file')
        [fail, result] = system(['export FSLDIR=' fsldir '; export PATH=' fsldir '/bin/:$PATH; source ' fsldir '/etc/fslconf/fsl.sh; ' fsldir '/bin/mcflirt -in ' data ' -refvol 0 -plots']);
        if fail
            disp(['WARNING: could not motion correct ' data])
        end
    end
    
    thisrun_params = load(paramsfile);
    thisrun_rot = thisrun_params(:,1:3);
    thisrun_rot_mm = thisrun_rot * brainradius;
    thisrun_params(:,1:3) = thisrun_rot_mm;
    thisrun_params_delta = [zeros(1,6) ; [thisrun_params(2:end,:) - thisrun_params(1:end-1,:)]];
    thisrun_params_delta_sq = thisrun_params_delta .^2;
    thisrun_FD = sum([sqrt(sum(thisrun_params_delta_sq(:,1:3),2)) sqrt(sum(thisrun_params_delta_sq(:,4:6),2))],2);
    
    FD = [FD ; thisrun_FD];
    
    if any(strcmp(thisseq_files(f).name,newfiles))
        newFD = [newFD ; thisrun_FD];
    end
    
    delete([data(1:end-7) '_mcf.nii.gz'])
    
end



tmask = (newFD <= FDthresh);

chunkedges = [1; diff([tmask;0])];
chunkstarts = find(chunkedges==1); chunkstops = (find(chunkedges==-1) -1);
chunksizes = zeros(length(chunkstarts),1);
chunkID = zeros(size(tmask));
for chunknum = 1:length(chunkstarts)
    chunkID(chunkstarts(chunknum) : chunkstops(chunknum)) = chunknum;
    chunksizes(chunknum) = nnz(chunkID==chunknum);
    if chunksizes(chunknum) < contiguoustimepoints
        tmask(chunkID==chunknum) = 0;
    end
end

pct_retained = nnz(tmask) / numel(tmask);
disp(['New runs: ' num2str(pct_retained*100) '% data retained from ' num2str(numel(tmask)) ' timepoints'])

plotted_tmask = single(tmask);
plotted_tmask(tmask==0) = NaN;

indices = [1:length(newFD)]';
figure; plot(indices,newFD,'-r',indices,repmat(FDthresh,length(newFD),1),'--k',indices,(plotted_tmask-1),'.g')
title('New runs')




if length(FD) > length(newFD)
    
    tmask = (FD <= FDthresh);
    
    chunkedges = [1; diff([tmask;0])];
    chunkstarts = find(chunkedges==1); chunkstops = (find(chunkedges==-1) -1);
    chunksizes = zeros(length(chunkstarts),1);
    chunkID = zeros(size(tmask));
    for chunknum = 1:length(chunkstarts)
        chunkID(chunkstarts(chunknum) : chunkstops(chunknum)) = chunknum;
        chunksizes(chunknum) = nnz(chunkID==chunknum);
        if chunksizes(chunknum) < contiguoustimepoints
            tmask(chunkID==chunknum) = 0;
        end
    end
    
    pct_retained = nnz(tmask) / numel(tmask);
    disp(['All runs: ' num2str(pct_retained*100) '% data retained from ' num2str(numel(tmask)) ' timepoints'])
    
    plotted_tmask = single(tmask);
    plotted_tmask(tmask==0) = NaN;
    
    indices = 1:length(FD);
    figure; plot(indices,FD,'-r',indices,repmat(FDthresh,length(FD),1),'--k',indices,(plotted_tmask-1),'.g')
    title('All runs')
    
end







