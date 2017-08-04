function fcp_cohortselector(datafile,targetdir,vcidfile,varargin)
%
% Name:fcp_cohortselector.m
% $Revision:$
% $Date:$
%
% jdp 3/2/12
%
% fcp_cohortselector(datafile,targetdir,vcidfile,*tmasklist*,**trimsize**)
%
% This script takes data that has undergone fcp_initial_process and helps
% the user sift through it to form a cohort of interest. 
%
% The datalist and targetdir are the same used in fcp_initial_process. The
% vcidfile is the same file used in the previous version of these scripts,
% which specifies the sex, age, twin status, and unique ID of subjects.
% This file has a single header line, and has 5 columns per subject:
% vcnum age sex uID twinID
% e.g.:
% vc11111 7.87  M 1 1
% vc22222 11.21 F 2 2
% vc33333 13.70 F 3 3
% vc44444 11.32 F 2 2
% vc55555 16.00 M 4 2
% 
%
% For this file, vc22222 was rescanned as subject vc44444 (hence they are
% both UID 2), and has a twin in vc55555 (hence the 2 in the twinUID for
% this separate person). The script will use this information to prevent
% twins or multiple vcnumbers representing a single person from entering
% analyses (if the user so chooses).
%  
% If a filename is passed in instead of a targetdir, the script reads the
% file to obtain specific targetdirs for each subject (if you are using a
% list of subjects whose fc data are stored in multiple places). 
% e.g.:
% vc11111 /data/place1/
% vc22222 /data/place1/
% vc33333 /data/place2/
%
% In general, the script will present a series of parameters and options to
% modify the parameters. Once parameters have been modified, the user can
% apply the settings to the datalist to see which subjects qualify under
% the settings. If the user likes the results, the cohort information can
% be printed out as several files, specified below.
% 
% OPTIONS:
% 11: FDthreshold - mask out volumes with FD>threshold
%     FDforward - mask out X volumes after a flagged volume
%     FDbackward - mask out X volumes before a flagged volume
% 12: DVthreshold - mask out volumes with DV>threshold
%     DVforward - mask out X volumes after a flagged volume
%     DVbackward - mask out X volumes before a flagged volume
%     DVtype - 1:use DV_333 (from fMRI processed data)
%              2:use DV_final (from fcprocessed data)
% 13: maskcombination - combine FD and DV masks with 1:AND or 2:OR
% 2:  agelo - minimum subject age for inclusion
%     agehi - maximum subject age for inclusion
% 3:  frameslo - minimum number of frames for inclusion
%     frameshi - maximum number of frames for inclusion
% 41: twins: 'same' considers twins as the same person
%            'separate' considers them to be different people
% 42: repeats: 'all' takes all scans from a person
%              'singlebest' takes the lowest DV scan from a person
%              'only' takes only subjects with repeat scans
% 77: apply the parameters to the cohort
%     This causes the following data for qualifying subject to be shown
%       M/F
%       Age
%       Volumes
%       % data lost
%       FD
%       DV
%       MVM
%       ddtMVM
%     A picture of the volumes, % data lost, FD, and DV is also made with
%     age on the X-axis and males in blue and females in red.
% 88: write out files for the cohort
%     This prompts for an outputname and writes the following files:
%       output/output_subjectproperties.txt
%       output/output_subjectproperties.tiff
%       output/output_datalist.txt (includes only qualifying subjects)
%       output/output_tmasklist.txt
%       output/tmask/vcXXXXX_output_tmask.txt
% 99: quit this script
%
% After running this script, you can run fcp_fcprocess.m
% 
% Suggested settings (superscrub) are:
%   FD<0.3, forward 3, backward 3
%   DVtype = DV_final
%   DV<3, forward 0, backward 0
%   OR operation
%
% Say you modify temporal masks somehow, or have some externally derived
% set of masks. If you wish to know the precise properties of the cohort
% paired with these masks, this script can operate on a supplied set of
% masks to generate the QC properties of the cohort.
% 
% If an optional tmasklist is supplied (one of the outputs of this script,
% listing a set of temporal masks), the script will bypass all of the
% options and calculate the QC properties of the cohort based upon the
% temporal masks in the tmasklist. The tmasklist is a txt file with 2
% columns, listing vcnumbers and hard paths to temporal masks, e.g.:
% vc11111 /data/src/tmask/vc11111_tmask.txt
% vc22222 /data/src/tmask/vc22222_tmask.txt
% vc33333 /data/src/tmask/vc33333_tmask.txt
% 
% This optional use will cause the output of the files:
%       output/output_subjectproperties.txt
%       output/output_subjectproperties.tiff
%
% If an additional optional integer is passed in, the script will ask the
% user for additional input, and will trim the tmasks in the tmasklist down
% to the trimsize, using either 1) the first X volumes, 2) a random set of 
% X volumes, or 3) the best X volumes by some QC measure.
% 
% fcp_cohortselector(datafile,targetdir,vcidfile,*tmasklist*,**trimsize**)
% fcp_cohortselector('data.txt','/data/is/here/','vcid.txt')
% fcp_cohortselector('data.txt','datalocations.txt','vcid.txt')
% fcp_cohortselector('data.txt','/data/is/here/','vcid.txt','tmasks.txt')
% fcp_cohortselector('data.txt','/data/is/here/','vcid.txt','tmasks.txt',250)
%
% NOTES:
% 5/24/12: added capability of multiple targetdirs
%          added uniqueIDs to outputs of cohort properties

% read in the subject data including location, vcnum, and boldruns
[dfile.prepdir dfile.prepstem dfile.paramname dfile.TR dfile.TRskip] = textread(datafile,'%s%s%s%f%d');
numsubs=size(dfile.prepdir,1);

if exist(targetdir)==7 % if it is a directory
    targetdir=repmat({targetdir},[numsubs 1]);
else
    [targetvcnum targetdir]=textread(targetdir,'%s%s');
    if ~isequal(dfile.prepstem,targetvcnum)
        error('Targetfile doesn''t match the datalist');
    end
end

if ~isempty(varargin)
	[vcs tmaskfile] = textread(varargin{1,1},'%s%s');
	if ~isequal(vcs,dfile.prepstem)
		error('VC ordering in datafile and tmasklist are not the same.');
	end
end

% read in the pertinent information
for i=1:numsubs
    subqc(i).srcdir=[targetdir{i,1} '/' dfile.prepstem{i,1} '/' ];
    subqc(i).FD=load([subqc(i).srcdir 'total_FD.txt']);
    subqc(i).DV_333=load([subqc(i).srcdir 'total_DV_333.txt']);
    subqc(i).DV_final=load([subqc(i).srcdir 'total_DV_final.txt']);
    subqc(i).mvm=load([subqc(i).srcdir 'total_movement_detrended.txt']);
    subqc(i).ddt_mvm=load([subqc(i).srcdir 'total_ddt_movement_detrended.txt']);
    subqc(i).orig_tmask=load([subqc(i).srcdir 'total_tmask.txt']);
    subqc(i).runborders=load([subqc(i).srcdir 'runborders.txt']);
    subqc(i).skipmask=ones(subqc(i).runborders(end,end),1);
    for k=1:size(subqc(i).runborders,1)
        subqc(i).skipmask(subqc(i).runborders(k,2):subqc(i).runborders(k,2)+dfile.TRskip(i)-1)=0;
    end
end

% read in the vcid file
[vcid.vc vcid.age vcid.sex vcid.uniqueID vcid.twinstogetheruniqueID] = vcIDfilereader(dfile.prepstem,vcidfile);

finished=0;
[params] = resetparams(vcid,subqc); 

if isempty(varargin) 
    
    while ~finished
        displaysettings(params);
        offeroptions();
        opt=input('Option: ');
        fprintf('\n');
        switch opt
            case 0
                [params] = resetparams(vcid,subqc);
            case 11
                [params.FDthresh]=input('Enter an FD threshold: ');
                [params.FDforward]=input('Enter TRs forward to augment: ');
                [params.FDbackward]=input('Enter TRs backward to augment: ');
            case 12
                badinput=1;
                while badinput
                    [params.DVtype]=input('Use DV values from pre-fc 333 (1) or final fc image (2): ');
                    switch params.DVtype
                        case {1,2}
                            badinput=0;
                    end
                end
                [params.DVthresh]=input('Enter a DV threshold: ');
                [params.DVforward]=input('Enter TRs forward to augment: ');
                [params.DVbackward]=input('Enter TRs backward to augment: ');
            case 13
                badinput=1;
                while badinput
                    [params.maskcombination]=input('Use AND (1) or OR (2) to combine FD and DV masks: ');
                    switch params.DVtype
                        case {1,2}
                            badinput=0;
                    end
                end
            case 14
                [params.snipsz]=input('What minimum number of contiguous TRs needed? (5 recommended): ');
                [params.runminsize]=input('What minimum number of volumes in run needed? (50 recommended): ');
            case 2
                [params.agelo]=input('Enter the minimum subject age desired: ');
                [params.agehi]=input('Enter the maximum subject age desired: ');
            case 3
                [params.frameslo]=input('Enter the lowest number of frames you wish to include: ');
                [params.frameshi]=input('Enter the highest number of frames you wish to include: ');
            case 41
                badinput=1;
                while badinput
                    [params.twins]=input('Do you want twins considered the same person or separate? (same/separate): ','s');
                    if isequal(params.twins,'same') || isequal(params.twins,'separate')
                        badinput=0;
                    else
                        fprintf('Please enter ''same'' or ''separate''\n');
                    end
                end
            case 42
                badinput=1;
                while badinput
                    [params.repeats]=input('Do you want all repeat scans, no repeat scans (take lowest FD RMS), or only repeat scans? (all/singlebest/only): ','s');
                    if isequal(params.repeats,'all') || isequal(params.repeats,'singlebest') || isequal(params.repeats,'only')
                        badinput=0;
                    else
                        fprintf('Please enter ''all'' or ''singlebest'' or ''only''\n');
                    end
                end
            case 77
                [verdict] = applyparams(params,subqc,vcid,dfile,0);
            case 88
                [verdict] = applyparams(params,subqc,vcid,dfile,1);
            case 99
                finished=1;
            otherwise
                fprintf('Try harder to follow instructions.\n');
                
        end
    end
    
else
    if numel(varargin)==1
        tmasklist=varargin{1,1};
        applyparams(params,subqc,vcid,dfile,1,tmasklist);
    elseif numel(varargin)==2
        tmasklist=varargin{1,1};
        trimfr=varargin{1,2};
        applyparams(params,subqc,vcid,dfile,1,tmasklist,trimfr);
    end
    
end
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This just shows the user their options

function offeroptions()

fprintf('\nChoose a parameter to modify:\n');
fprintf('0:  Reset all parameters\n');
fprintf('11: FD threshold\n');
fprintf('12: DV threshold\n');
fprintf('13: Mask AND/OR\n');
fprintf('14: Minimum Frames specification:\n');
fprintf('2:  Age\n');
fprintf('3:  Frames\n');
fprintf('41: Twins\n');
fprintf('42: Repeats\n');
% fprintf('7:  Store properties and qualifying subjects for comparisons\n');
% fprintf('8:  Perform comparisons\n');
fprintf('77: Apply settings\n');
fprintf('88: Print files for settings\n');
fprintf('99: Finished\n');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params] = resetparams(vcid,subqc)

params.FDthresh=max(arrayfun(@(x)(max(x.FD)),subqc));
params.FDforward=0;
params.FDbackward=0;
params.DVtype=2; % post-fc is default
params.DVthresh=max(arrayfun(@(x)(max(x.DV_final)),subqc));
params.DVforward=0;
params.DVbackward=0;
params.maskcombination=1; % AND is default
params.runminsize = 50;
params.snipsz = 5;
params.FRAMElimit=max(arrayfun(@(x)(numel(x.DV_final)),subqc));
params.agelo=min(vcid.age);
params.agehi=max(vcid.age);
params.frameslo=0;
params.frameshi=max(arrayfun(@(x)(numel(x.DV_final)),subqc));
params.twins='same'; % treat twins as identical people
params.repeats='singlebest'; % do not allow people scanned under multiple vcnums


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function displaysettings(params)

fprintf('\n\n####SETTINGS####\n');
fprintf('FDthreshold\t%g\n',params.FDthresh);
fprintf('FDforward\t%d\n',params.FDforward);
fprintf('FDbackward\t%d\n',params.FDbackward);
fprintf('DVtype\t\t%d (1=pre-fc, 2=post-fc)\n',params.DVtype);
fprintf('DVthreshold\t%g\n',params.DVthresh);
fprintf('DVforward\t%d\n',params.DVforward);
fprintf('DVbackward\t%d\n',params.DVbackward);
fprintf('maskcombination\t%d (1=AND, 2=OR)\n',params.maskcombination);
fprintf('FRAMElimit\t%d\n',params.FRAMElimit);
fprintf('Agelo\t\t%g\n',params.agelo);
fprintf('Agehi\t\t%g\n',params.agehi);
fprintf('Frameslo\t%d\n',params.frameslo);
fprintf('Frameshi\t%d\n',params.frameshi);
fprintf('Twins\t\t%s (same treats twins as the same person)\n',params.twins);
fprintf('Repeats\t\t%s (singlebest takes the best scan of a person scanned under multiple vcnumbers)\n',params.repeats);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [verdict] = applyparams(params,subqc,vcid,dfile,writeoutput,varargin)

externaltmask=0;
trimmasks=0;
if ~isempty(varargin)
    externaltmask=1;
end
if nargin==7
    trimmasks=1;
end


if ~externaltmask % if generating tmasks from scratch
    
    % form the FD masks
    tempfd=arrayfun(@(x)(x.FD<=params.FDthresh),subqc,'uniformoutput',0);
    tempfdaug=cellfun(@(x)fcp_augmenttmask(x,params.FDbackward,params.FDforward),tempfd,'uniformoutput',0);
    
    % form the DV masks
    if params.DVtype==1
        tempdv=arrayfun(@(x)(x.DV_333<=params.DVthresh),subqc,'uniformoutput',0);
    else
        tempdv=arrayfun(@(x)(x.DV_final<=params.DVthresh),subqc,'uniformoutput',0);
    end
    tempdvaug=cellfun(@(x)fcp_augmenttmask(x,params.DVbackward,params.DVforward),tempdv,'uniformoutput',0);
    
    % combine the masks
    if params.maskcombination==1
        combo=cellfun(@(x,y)~(~x&~y),tempfdaug,tempdvaug,'uniformoutput',0);
    else
        combo=cellfun(@(x,y)~(~x|~y),tempfdaug,tempdvaug,'uniformoutput',0);
    end
    
    % apply the skipmask
    junkst=(cell2struct(combo,'combo',1))';
    combo=arrayfun(@(x,y)x.combo&y.skipmask,junkst,subqc,'uniformoutput',0);
    
    % apply minimum number of frames and continuous frames mask 
    for i=1:size(subqc,2) % cycle subjects
        combo{1, i} = sniptinymask(combo{1, i}, subqc(i).runborders ,params.snipsz);  
    end
    
    
    % THERE SHOULD BE A CHECK/FIX FOR EMPTY RUNS OR RUNS WITH <X TRS!!!!!!
    %runminsize=50;
    runminsize = params.runminsize;
    for i=1:size(subqc,2) % cycle subjects
    	finalcombo{1,i}=zeros(numel(combo{1,i}),1);
        for k=1:size(subqc(i).runborders,1) % cycle runs
            r=nnz(combo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3)));
            if r<runminsize
                combo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3))=0;
                fprintf('8====> %s run%d has %d volumes; need %d; excluding (:-(\n',vcid.vc{i},k,r,runminsize);
                keptrun{i}(k,1)=0;
                finalcombo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3))=0;
            else
                keptrun{i}(k,1)=1;
                finalcombo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3))=1;
            end
        end
    end
    
else % if tmasks are provided by the user
    [vcs tmaskfile] = textread(varargin{1,1},'%s%s');
    for i=1:numel(vcs)
        combo{1,i}=logical(load(tmaskfile{i,1}));
    end
    
    if trimmasks % if they want to trim masks to identical numbers of frames
        %runminsize=50;
        runminsize = params.runminsize;
        combo=trimframes(combo,subqc,varargin{1,2},runminsize);
        % THERE SHOULD BE A CHECK/FIX FOR EMPTY RUNS OR RUNS WITH <X TRS!!!!!!
        
        for i=1:size(subqc,2) % cycle subjects
            finalcombo{1,i}=zeros(numel(combo{1,i}),1);
            for k=1:size(subqc(i).runborders,1) % cycle runs
                r=nnz(combo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3)));
                if r<runminsize
                    combo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3))=0;
                    fprintf('8====> %s run%d has %d volumes; need %d; excluding (:-(\n',vcid.vc{i},k,r,runminsize);
                    keptrun{i}(k,1)=0;
                    finalcombo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3))=0;
                else
                    keptrun{i}(k,1)=1;
                    finalcombo{1,i}(subqc(i).runborders(k,2):subqc(i).runborders(k,3))=1;
                end
            end
        end
    end
end

% calculate within-mask FD and DV for repeat/twin considerations
st=(cell2struct(combo,'combo',1))';
maskFD=arrayfun(@(x,y)array_calc_rms(x.FD(y.combo)),subqc,st);
maskDV_333=arrayfun(@(x,y)array_calc_rms(x.DV_333(y.combo)),subqc,st);
maskDV_final=arrayfun(@(x,y)array_calc_rms(x.DV_final(y.combo)),subqc,st);
if params.DVtype==1
    maskDV=maskDV_333;
else
    maskDV=maskDV_final;
end

% calculate within-mask movement parameters
maskMVM=arrayfun(@(x,y)array_calc_rms(x.mvm(y.combo,:)),subqc,st);
maskddtMVM=arrayfun(@(x,y)array_calc_rms(x.ddt_mvm(y.combo,:)),subqc,st);

% decide whether sufficient volumes remain
frames.all=cellfun(@numel,combo);
frames.kept=cellfun(@nnz,combo);
frames.lost=frames.all-frames.kept;
frames.propkept=frames.kept./frames.all;

% apply the twins and repeat exclusion criteria
verdict.enoughframes=frames.kept'>=params.frameslo & frames.kept'<=params.frameshi;
verdict.ageisgood=vcid.age>=params.agelo & vcid.age<=params.agehi;
verdict.ageandframes = verdict.enoughframes & verdict.ageisgood;
verdict.finalmask = logical(twinrepeat_remover(verdict.ageandframes,maskDV,params.twins,params.repeats,vcid.uniqueID,vcid.twinstogetheruniqueID));
% if tmasks already provided ensure that no subjects are actually lost
if ~isempty(varargin)
    verdict.finalmask=verdict.finalmask | ~verdict.finalmask; % set all values to 1
end

% determine cohort properties
[tempmale tempfemale sexrgb] = sexcounter(vcid.sex(verdict.finalmask));
[tempagemin tempagemax tempagemean tempagestd] = minmaxmean(vcid.age(verdict.finalmask));
[tempframesallmin tempframesallmax tempframesallmean tempframesallstd] = minmaxmean(frames.all(verdict.finalmask));
[tempframeslostmin tempframeslostmax tempframeslostmean tempframesloststd] = minmaxmean(frames.lost(verdict.finalmask));
[tempframeskeptmin tempframeskeptmax tempframeskeptmean tempframeskeptstd] = minmaxmean(frames.kept(verdict.finalmask));
[temppropframesmin temppropframesmax temppropframesmean temppropframesstd] = minmaxmean(frames.propkept(verdict.finalmask));
[tempFDmin tempFDmax tempFDmean tempFDstd] = minmaxmean(maskFD(verdict.finalmask));
[tempDVmin tempDVmax tempDVmean tempDVstd] = minmaxmean(maskDV(verdict.finalmask));
[tempDV333min tempDV333max tempDV333mean tempDV333std] = minmaxmean(maskDV_333(verdict.finalmask));
[tempDVfinalmin tempDVfinalmax tempDVfinalmean tempDVfinalstd] = minmaxmean(maskDV_final(verdict.finalmask));
[tempMVMmin tempMVMmax tempMVMmean tempMVMstd] = minmaxmean(maskMVM(verdict.finalmask));
[tempddtMVMmin tempddtMVMmax tempddtMVMmean tempddtMVMstd] = minmaxmean(maskMVM(verdict.finalmask));

%%% display the results
if writeoutput==0
    fprintf('\n\n####COHORT PROPERTIES####\n')
    fprintf('Qualifying subjects:\t%d/%d\n',nnz(verdict.finalmask),numel(verdict.finalmask));
    fprintf('Male/Female:\t\t%d/%d\n',tempmale,tempfemale);
    fprintf('Age:\t\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempagemean,tempagestd,tempagemin,tempagemax);
    fprintf('Volumes:\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempframeskeptmean,tempframeskeptstd,tempframeskeptmin,tempframeskeptmax);
    fprintf('%% data lost:\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',temppropframesmean,temppropframesstd,temppropframesmin,temppropframesmax);
    fprintf('RMS FD:\t\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempFDmean,tempFDstd,tempFDmin,tempFDmax);
    fprintf('RMS DV (chosen):\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempDVmean,tempDVstd,tempDVmin,tempDVmax);
    fprintf('RMS DV_333:\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempDV333mean,tempDV333std,tempDV333min,tempDV333max);
    fprintf('RMS DV_final:\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempDVfinalmean,tempDVfinalstd,tempDVfinalmin,tempDVfinalmax);
    fprintf('RMS MVM:\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempMVMmean,tempMVMstd,tempMVMmin,tempMVMmax);
    fprintf('RMS ddtMVM:\t\t%4.2f(%4.2f)\t(%4.2f-%4.2f)\n',tempddtMVMmean,tempddtMVMstd,tempddtMVMmin,tempddtMVMmax);
    if nnz(verdict.finalmask) > 0
        plotrms(sexrgb,vcid.age(verdict.finalmask),frames.kept(verdict.finalmask),frames.propkept(verdict.finalmask),maskFD(verdict.finalmask),maskDV(verdict.finalmask));
    end
    fprintf('Hit a key');
    pause;
    close;
    
else
    
    startdir=pwd;
    outname=input('Enter a name to tag output files: ','s');
    if ~exist(outname, 'dir')
        mkdir(outname);
    end
    
    % should output be sorted?
    positions=find(verdict.finalmask);
    sortby=0;
    if ~externaltmask
        sortby=input('Do you wish to order output? (0=no,1=age,2=frames,3=propframesleft,4=FD,5=DV,6=MVM) ');
    end
    if sortby
        notfinished=1;
        while notfinished
            switch sortby
                case 1
                    sortmat=vcid.age(verdict.finalmask);
                    notfinished=0;
                case 2
                    sortmat=[frames.kept(verdict.finalmask)]';
                    notfinished=0;
                case 3
                    sortmat=[frames.propkept(verdict.finalmask)]';
                    notfinished=0;
                case 4
                    sortmat=[maskFD(verdict.finalmask)]';
                    notfinished=0;
                case 5
                    sortmat=[maskDV(verdict.finalmask)]';
                    notfinished=0;
                case 6
                    sortmat=[maskMVM(verdict.finalmask)]';
                    notfinished=0;
                otherwise
            end
        end
        
        order=input('Sort (1) ascending or (2) descending? ');
        notfinished=1;
        while notfinished
            switch order
                case 1
                    sorttype='ascend';
                    notfinished=0;
                case 2
                    sorttype='descend';
                    notfinished=0;
                otherwise
                    error('Have to use 1 or 2 for sorting.');
            end
        end
        [s ind] = sort(sortmat,1,sorttype);
    else
        ind=[1:nnz(verdict.finalmask)]';
    end
    outputorder=positions(ind);
    
    
    % write data for each subject
    fid=fopen([outname '/' outname '_subjectproperties.txt'],'w');
    fprintf(fid,'VC\tUID\ttwinUID\tSex\tAge\tFrames\tFkept\tFlost\t%%FRkept\tRMS_FD\tRMS_DV\tRMS_DV_333\tRMS_DV_final\tRMS_MVM\tRMS_ddtMVM\n');
    for i=1:numel(outputorder)
        fprintf(fid,'%s\t%d\t%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n',vcid.vc{outputorder(i)},vcid.uniqueID(outputorder(i)),vcid.twinstogetheruniqueID(outputorder(i)),vcid.sex{outputorder(i)},vcid.age(outputorder(i)),frames.all(outputorder(i)),frames.kept(outputorder(i)),frames.lost(outputorder(i)),frames.propkept(outputorder(i)),maskFD(outputorder(i)),maskDV(outputorder(i)),maskDV_333(outputorder(i)),maskDV_final(outputorder(i)),maskMVM(outputorder(i)),maskddtMVM(outputorder(i)));
    end
    fprintf(fid,'\n');
    fprintf(fid,'\t%dM\t%4.2g\t%5.3g\t%5.3g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\n',tempmale,tempagemean,tempframesallmean,tempframeskeptmean,tempframeslostmean,temppropframesmean,tempFDmean,tempDVmean,tempDV333mean,tempDVfinalmean,tempMVMmean,tempddtMVMmean);
    fprintf(fid,'\t%dF\t%4.2g\t%5.3g\t%5.3g\t%4.2g\t%g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\n',tempfemale,tempagestd,tempframesallstd,tempframeskeptstd,tempframesloststd,temppropframesstd,tempFDstd,tempDVstd,tempDV333std,tempDVfinalstd,tempMVMstd,tempddtMVMstd);
    fprintf(fid,'\t\t%g\t%g\t%5.3g\t%5.3g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\n',tempagemin,tempframesallmin,tempframeskeptmin,tempframeslostmin,temppropframesmin,tempFDmin,tempDVmin,tempDV333min,tempDVfinalmin,tempMVMmin,tempddtMVMmin);
    fprintf(fid,'\t\t%g\t%g\t%5.3g\t%5.3g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\t%4.2g\n',tempagemax,tempframesallmax,tempframeskeptmax,tempframeslostmax,temppropframesmax,tempFDmax,tempDVmax,tempDV333max,tempDVfinalmax,tempMVMmax,tempddtMVMmax); 
    if ~externaltmask
        fprintf(fid,'\n\n\n\n');
        fprintf(fid,'FDthresh:\t%g\n',params.FDthresh);
        fprintf(fid,'FDforward:\t%d\n',params.FDforward);
        fprintf(fid,'FDbackward:\t%d\n',params.FDbackward);
        fprintf(fid,'DVtype:\t\t%d\t\t(1=DV_333,2=DV_final)\n',params.DVtype);
        fprintf(fid,'DVthresh:\t%g\n',params.DVthresh);
        fprintf(fid,'DVforward:\t%d\n',params.DVforward);
        fprintf(fid,'DVbackward:\t%d\n',params.DVbackward);
        fprintf(fid,'maskcomb:\t%d\t\t(1=AND,2=OR)\n',params.maskcombination);
        fprintf(fid,'agelo:\t\t%g\n',params.agelo);
        fprintf(fid,'agehi:\t\t%g\n',params.agehi);
        fprintf(fid,'frameslo:\t%g\n',params.frameslo);
        fprintf(fid,'frameshi:\t%g\n',params.frameshi);
        fprintf(fid,'twins:\t\t%s\t\t(''same'' treats twins as one person)\n',params.twins);
        fprintf(fid,'repeats:\t%s\t(''singlebest'' takes the lowest DV session of a person)\n',params.repeats);
    end
    fclose(fid);
    if nnz(outputorder)>0
        plotrms(sexrgb,vcid.age(verdict.finalmask),frames.kept(verdict.finalmask),frames.propkept(verdict.finalmask),maskFD(verdict.finalmask),maskDV(verdict.finalmask));
        saveas(gcf,[outname '/' outname '_subjectproperties.tiff'],'tiff');
    end
    close;
    
   
    if (~externaltmask | trimmasks)
        % write the temporal masks
        if ~exist([startdir '/' outname '/tmask/'])
            mkdir([startdir '/' outname '/tmask/']);
        end
        for i=1:numel(outputorder)
        	if ((~externaltmask) | trimmasks)
         	   		reduced_tmaskname{i,1}=[startdir '/' outname '/tmask/' outname '_' vcid.vc{outputorder(i)} '_reduced_tmask.txt'];
         	   		dlmwrite(reduced_tmaskname{i,1},combo{outputorder(i)}(logical(finalcombo{1,outputorder(i)})));
                    tmaskname{i,1}=[startdir '/' outname '/tmask/' outname '_' vcid.vc{outputorder(i)} '_tmask.txt'];
                    dlmwrite(tmaskname{i,1},combo{outputorder(i)});
            end

            
            % if need to print out modified paramfiles
            if nnz(~keptrun{1,outputorder(i)})>0
                newprms=subqc(outputorder(i)).runborders(logical(keptrun{1,outputorder(i)}),1);
                dfile.paramname{outputorder(i)}=[startdir '/' outname '/tmask/' vcid.vc{outputorder(i)} '_reduced.params' ];
                fid=fopen(dfile.paramname{outputorder(i)},'w');
                fprintf(fid,'set boldruns = (');
                for k=1:numel(newprms)
                    fprintf(fid,'%d',newprms(k));
                    if k~=numel(newprms)
                        fprintf(fid,'\t');
                    else
                        fprintf(fid,')\n');
                    end
                end
            end
                
        end
        
        % write the reduced temporal mask list
        fid=fopen([outname '/' outname '_reduced_tmasklist.txt'],'w');
        for i=1:numel(outputorder)
            fprintf(fid,'%s\t%s\n',vcid.vc{outputorder(i)},reduced_tmaskname{i,1});
        end
        fclose(fid);
        
        % write the non-reduced temporal mask list
        fid=fopen([outname '/' outname '_nonreduced_tmasklist.txt'],'w');
        for i=1:numel(outputorder)
            fprintf(fid,'%s\t%s\n',vcid.vc{outputorder(i)},tmaskname{i,1});
        end
        fclose(fid);
        
        % write out an updated datalist with only the interesting subjects
        fid=fopen([outname '/' outname '_reduced_datalist.txt'],'w');
        for i=1:numel(outputorder)
            fprintf(fid,'%s\t%s\t%s\t%2.2f\t%d\n',dfile.prepdir{outputorder(i)},dfile.prepstem{outputorder(i)},dfile.paramname{outputorder(i)},dfile.TR(outputorder(i)),dfile.TRskip(outputorder(i)));
        end
        fclose(fid);
    end
    
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotrms(sexrgb,age,frames,propframeskept,FD,DV)

subplot(2,2,1);
hold on;
for i=1:numel(age)
    plot(age(i),frames(i),'.','Color',sexrgb(i,:));
end
hold off;
ylabel('Frames'); xlabel('Age'); ylim([0 max(frames)+100]);
subplot(2,2,2);
hold on;
for i=1:numel(age)
    plot(age(i),propframeskept(i),'.','Color',sexrgb(i,:));
end
hold off;
ylabel('% data left'); xlabel('Age'); ylim([0 1]);
subplot(2,2,3);
hold on;
for i=1:numel(age)
    plot(age(i),FD(i),'.','Color',sexrgb(i,:));
end
hold off;
ylabel('RMS FD'); xlabel('Age'); ylim([0 max(FD)+.5]);
subplot(2,2,4);
hold on;
for i=1:numel(age)
    plot(age(i),DV(i),'.','Color',sexrgb(i,:));
end
hold off;
ylabel('RMS DV'); xlabel('Age');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This function applies the repeats and twins settings to a cohort of
% subjects who otherwise meet the parameter space. A winnowed cohort of
% qualifying subjects is returned.

function [parametermask] = twinrepeat_remover(parametermask,temprms,twins,repeats,uniqueID,twinstogetheruniqueID)

% are twins considered the same or separate persons?
switch twins
    case 'same'
        ID=twinstogetheruniqueID;
    case 'separate'
        ID=uniqueID;
end

% who are the unique people in the qualifying ID list?

tempID=ID.*parametermask; % mask out subjects who don't qualify
uniques=unique(tempID); % find the unique ID values
if uniques(1)==0 % if subjects have been masked out, ignore their tempID value (0)
    uniques(1)=[];
end
d=size(uniques);

if ~isempty(uniques) % it is possible nobody qualifies anymore, don't make this crash the program
    
    % now winnow according to the repeats settings
    switch repeats
        case 'all' % do nothing
        case 'only' % the ID must appear >=2 times to qualify
            for i=1:d(1)
                if nnz(ismember(tempID,uniques(i)))==1
                    parametermask(ismember(tempID,uniques(i)))=0;
                end
            end
        case 'singlebest' % only the lowest RMS scan is retained from a person, other scans are discarded
            for i=1:d(1)
                if nnz(ismember(tempID,uniques(i)))>1
                    overlaps=ismember(tempID,uniques(i));
                    maskedrms=overlaps.*temprms';
                    maskedrms(maskedrms==0)=NaN;
                    [minrms minposition]=min(maskedrms);
                    overlaps(minposition)=0;
                    parametermask=parametermask-overlaps;
                end
            end
        otherwise 
            fprintf('Please enter ''all'' ''only'' or ''singlebest'' \n');
    end
    
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This just calculates the minimum, maximum, mean, and standard deviation
% of some array that is passed in to the function.

function [somemin somemax somemean somestd] = minmaxmean(somemat)

somemin=min(somemat(:));
somemax=max(somemat(:));
somemean=mean(somemat(:));
somestd=std(somemat(:));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This calculates the number of males and females from a cell array of M/F

function [male female sexrgb] = sexcounter(sex)

sexrgb=[];
numsubs=size(sex,1);
male=0;
female=0;
for i=1:numsubs
    switch sex{i,1}
        case 'M'
            male=male+1;
            sexrgb=[sexrgb;0 0 1];
        case {'F','female'}
            female=female+1;
            sexrgb=[sexrgb;1 0 0];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This trims masks to a particular number of volumes

function [combo]=trimframes(combo,subqc,X,runmin)

numsubs=numel(combo);
for i=1:numsubs
    
    fprintf('Trimming %d %s\n',i,subqc(i).srcdir);
    [a]=find(combo{1,i});
    if (numel(a)<X)                                                         % if not enough good TRs
        fprintf('Subject %d doesn''t even have %d good frames\n',i,X);
    elseif (numel(a)==X)                                                    % if exactly enough TRs
    else
        
        % determine good TRs within each run
        b=subqc(i).runborders;  
        for j=1:size(b,1)
            b(j,4)=j; % index for resorting the data
            b(j,5)=nnz(a>=b(j,2)&a<=b(j,3)); % number of good TRs in the run
            b(j,6)=b(j,5)-runmin; % extra TRs that could be cut if needed
        end
        
        [b bindex]=sortrows(b,5); % sort by number of good TRs
        b=flipud(b); bindex=flipud(bindex);
        
        keepgoing=1;
        k=0;
        goodTRs=0;
        extraTRs=0;
        while keepgoing
            k=k+1;
            if k>size(b,1) % if index exceeds numruns
                keepgoing=0;
                break;
            end
            
            goodTRs=goodTRs+b(k,5);
            extraTRs=extraTRs+b(k,6);
            
            if goodTRs>=X % if sufficient TRs now identified
                if b(k,6)>=(goodTRs-X) & (runmin*k<=X) % if can just cut TRs in this most recent run
                    a=trimTRs(a,(goodTRs-X),b(k,2),b(k,3));
                    while k<size(b,1) % if not at the final run zero out others
                        k=k+1;
                        a=trimTRs(a,b(k,5),b(k,2),b(k,3));
                    end
                    if numel(a)~=X
                        error('SOMETHING WASN''T RIGHT WITH TRIMMING THE RUN');
                    end
                    keepgoing=0; %FINISHED
                    
                elseif (extraTRs>=(goodTRs-X)) & (runmin*k<=X) % if can cut TRs in this and previous runs to achieve X
                    kk=k;
                    while k<size(b,1) % if not at the final run zero out others
                        k=k+1;
                        a=trimTRs(a,b(k,5),b(k,2),b(k,3));
                    end
                    
                    for k=kk:-1:1
                        while (numel(a)>X) & (b(k,6)>0) % while needing to zero out TRs still to achieve X
                            a=trimTRs(a,1,b(k,2),b(k,3));
                            b(k,6)=b(k,6)-1;
                        end
                    end
                    
                    if numel(a)~=X
                        error('SOMETHING WASN''T RIGHT WITH TRIMMING THE RUN');
                    end
                    keepgoing=0; %FINISHED
                    
                elseif runmin*k>X
                    fprintf('Stuck with too many runs that can''t be trimmed to the right total length; person %d literally can''t be trimmed\n',i);
                    keepgoing=0;
                    
                % the only other possibility is that there weren't adequate
                % extraTRs to trim down to X AND minsize*k is not yet > X
                end
            end
        end
        
        % set the correct TRs to 1
        combo{1,i}=logical(combo{1,i}*0);
        combo{1,i}(a)=1;
        
        if nnz(combo{1,i})~=X
            error('SOMETHING WASN''T RIGHT WITH TRIMMING THE RUN');
        end
    end
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TRs]=trimTRs(TRs,numtotrim,startval,endval)

a=find(TRs>=startval & TRs<=endval);
if numel(a)>=numtotrim
    a(1:end-numtotrim)=[];
    TRs(a)=[];
else
    error('There are not enough values to trim in trimTRs');
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
function [mask2]=sniptinymask(mask1,borders,snipsz)

runs=size(borders,1);
for r=1:runs
    tempmask1=mask1(borders(r,2):borders(r,3));
    [chunksize startnums endnums] = maskgaps(~tempmask1);
    goodsize=endnums-startnums+1;
    goodlost=goodsize<=snipsz;
    for j=1:numel(goodlost)
        if goodlost(j)
            tempmask1(startnums(j):endnums(j))=0;
        end
    end
    mask2(borders(r,2):borders(r,3))=tempmask1;
    
end
    
mask2=mask2';
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [chunksize startnums endnums] = maskgaps(mask)

% presumes that mask is 1=keep, 0=discard


if nnz(mask)==0 % if all discarded
    chunksize=-numel(mask);
elseif nnz(mask)==numel(mask) % if all kept
    chunksize=numel(mask);
else % if some kept and some discarded
    
    dmask=diff(mask);
    numchunks=nnz(dmask)+1;
    chunksize=ones(numchunks,1);
    chunknum=1;
    
    for i=1:numel(dmask)
        if dmask(i)==0
            chunksize(chunknum,1)=chunksize(chunknum,1)+1;
        elseif dmask(i)==1
            chunksize(chunknum,1)=-chunksize(chunknum,1);
            chunknum=chunknum+1;
        elseif dmask(i)==-1
            chunknum=chunknum+1;
        end
        
        if i==numel(dmask)
            if chunksize(chunknum-1,1)>0
                chunksize(chunknum,1)=-chunksize(chunknum,1);
            end
        end
        
    end
end

if ~isequal((sum(abs(chunksize))),numel(mask))
    disp('error');
end

endnums=cumsum(abs(chunksize));
endnums=endnums(chunksize<0);
startnums=endnums+chunksize(chunksize<0)+1;
        
        
        
        
        