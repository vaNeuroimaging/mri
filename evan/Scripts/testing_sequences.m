sequences = {'RSFCSENSE','RSFCTR3SENSE'};%{'restingstate','RSFC_std','RSFC_TR3','RSFC_443','RSFC_444','RSFC_noSen_TR2p5','RSFC_noSen_444'};
subjects = {'Zach','Tanner','Joel'};%'Ramy',
brainradius = 50;
contiguoustimepoints = 5;
FDthresh = .2;

FDs = cell(1,0);
tmasks = cell(1,0);
meanFDs = [];
stdFDs = [];
pcts_retained = [];
names = cell(1,0);
seqs = [];

for s = 1:length(subjects)
    
    cd(['/home/data/subjects/' subjects{s} '/raw/'])
    
    for seq = 1:length(sequences)
        
%        filenames = dir([sequences{seq} '*.nii.gz']);
        filenames = dir([sequences{seq} '*_mcf.par']);
        
        for f = 1:length(filenames)
%             
%             system(['mcflirt -in ' filenames(f).name ' -refvol 0 -plots']);
%             
%             data = load_untouch_nii_2D([filenames(f).name(1:end-7) '_mcf.nii.gz']);
%             out = demean_detrend(data.img(:,5:end));
%             data.hdr.dime.dim(5) = data.hdr.dime.dim(5) - 4;
%             data.img = out;
%             save_untouch_nii_2D(data,[filenames(f).name(1:end-7) '_mcf_dmdt.nii.gz']);
%             system(['fslview ' filenames(f).name(1:end-7) '_mcf_dmdt.nii.gz -l render3 -b -15,15 &'])
%             
%            %delete([filenames(f).name(1:end-7) '_mcf.nii.gz']);
%            calc_FD([filenames(f).name(1:end-7) '_mcf.par'])
%            thisrun_params = load([filenames(f).name(1:end-7) '_mcf.par']);
            calc_FD([filenames(f).name])
            thisrun_params = load([filenames(f).name]);
            
            thisrun_rot = thisrun_params(:,1:3);
            thisrun_rot_mm = thisrun_rot * brainradius;
            thisrun_params(:,1:3) = thisrun_rot_mm;
            thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];
            
            %demean and detrend motion regressors
            thisrun_params_demeandetrend = demean_detrend(thisrun_params');
            thisrun_params_delta_demeandetrend = demean_detrend(thisrun_params_delta');
            
            %Calculate FD and tmask
            thisFD = [sum(abs(thisrun_params_delta),2)]';
            FDs(end+1) = {thisFD};
            
            meanFDs(end+1) = mean(thisFD);
            stdFDs(end+1) = std(thisFD);
            names(end+1) = {[subjects{s} ': ' sequences{seq} '_' num2str(f)]};
            
            
            thisrun_tmask = (thisFD') < FDthresh;
            
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
            
            tmasks(end+1) = {thisrun_tmask};
            pcts_retained(end+1) = nnz(thisrun_tmask) / numel(thisrun_tmask);
            seqs(end+1) = seq;
            
            
            
        end
        
    end
    
end
% figure
% errorbar(meanFDs,stdFDs,'k.')
% set(gca,'XTick',[1:length(meanFDs)])
% set(gca,'XTickLabel',names)
% set(gca,'Ylim',[min(meanFDs - stdFDs) - .01 , max(meanFDs + stdFDs) + .01])
% 
figure
b = bar(pcts_retained,'k');
set(gca,'XTick',[1:length(pcts_retained)])
set(gca,'XTickLabel',names)