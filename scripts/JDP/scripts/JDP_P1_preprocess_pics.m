function JDP_P1_preprocess_pics(varargin)

% This script helps check preprocessing
%
% Inputs: (pdir,scriptstr,startsub,endsub,analtype)
% pdir - the parent directory, a hard path
% scriptstr - the processing script name, determines the log directory
% startsub - first subject in cleansublist.txt to process
% endsub - final subject to process
% analtype - which analysis to perform
%   1: parses log files for errors and warnings
%   2: creates many .pngs for each subject (in log/rawimgs)
%   20: groups output of (2) to usefully examine particular steps

pdir=varargin{1,1};
scriptstr=varargin{1,2};
startsub=varargin{1,3};
if ischar(startsub)
    startsub=str2num(startsub);
end
endsub=varargin{1,4};
if ischar(endsub)
    endsub=str2num(endsub);
end
analtype=varargin{1,5};

% basics of the file structure
datadir=[pdir '/subs'];
scriptdir=[pdir '/scripts'];

% where we store summary logs and pictures
[jk logstem jk2]=fileparts(scriptstr);
logdir=['logs_' logstem ];
logpath=[scriptdir '/' logdir ];
system(['mkdir -p ' logpath]);

% the subjects that could be processed
subjlistfile=[scriptdir '/cleansublist.txt'];
subjlist=textread(subjlistfile,'%s');

switch analtype
    
    case 'assemblelogs'
        
        %%%
        %%% this parses log files to create summaries of errors/warnings
        %%%
        
        %this will concatenate all the subject warnings in a summary file
        ofile=[logpath '/SUMMARY_i' num2str(startsub) '_to_i' num2str(endsub) '_WARNING.txt'];
        system(['rm ' ofile]); system(['touch ' ofile]);
        for i=startsub:endsub
            sub=subjlist{i};
            for k=1:10
                str=['echo subject ' num2str(sub) ' >> ' ofile];
                system(str);
            end
            warnfile = [logpath '/logs/' num2str(sub) '_warning.txt'];
            str=['cat ' warnfile ' >> ' ofile];
            system(str);
        end
        
        % now pull only the errors out
        ofile=[logpath '/SUMMARY_i' num2str(startsub) '_to_i' num2str(endsub) '_ERROR.txt'];
        system(['rm ' ofile]); system(['touch ' ofile]);
        for i=startsub:endsub
            sub=subjlist{i};
            for k=1:1
                str=['echo subject ' num2str(sub) ' >> ' ofile];
                system(str);
            end
            warnfile = [logpath '/logs/' num2str(sub) '_warning.txt'];
            str=['cat ' warnfile ' | grep -n ''ERROR\|Error'' >> ' ofile ];
            system(str);
        end
        
        
    case 'basicpictures'
        
        % Pictures of preprocessing results (in /logdir/rawimgs)
        % Naming scheme:
        %   Subject position in cleansublist * 1Million
        %   0-10        MPRAGE skullstripping
        %   11-999      EPI skullstrippings
        %   1000-1009   MPR-ATL registration
        %   1010-1999   EPI-MPR registrations
        %   2000-2999   EPI-ATL registrations
        %   3000-3999   FS segmentations and erosions over MPRAGE
        %   4000-4999   FS segmentations and erosions over EPIs
        %   5000-5999   methods of EPI registration
        %   6000-6999   method of MPRAGE registration
        
        
        subsdir=[pdir '/subs'];
        
        file.Power264=[scriptdir '/Power264_stack.nii.gz'];
        clrs.Power264=jet(256);
        lims.Power264=[0 1];
        
        switch scriptstr
            case {'JDP_preprocess.sh','JDP_P1_preprocess.sh'}
                resdir='P1.results';
        end
        
        for i=startsub:endsub
            fprintf('Subject %d\n',i);
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            system(['rm -rf ' picpath ]);
            system(['mkdir -p ' picpath ]);
            
            
            % structural images of interest
            limz.tt_n27=[0 300];
            rgb.tt_n27=gray(256);
            file.TT_N27=[subdir '/TT_N27.nii.gz'];
            clrs.TT_N27=rgb.tt_n27;
            lims.TT_N27=limz.tt_n27;
            file.TT_N27_EPI=[subdir '/TT_N27_EPI.nii.gz'];
            clrs.TT_N27_EPI=rgb.tt_n27;
            lims.TT_N27_EPI=limz.tt_n27;
            
            % set scaling for the mprage
            file.mprage=[subdir '/mprage.nii.gz'];
            
            imthr=.9;
            tmp=load_untouch_nii(file.mprage);
            tmpimg=tmp.img(:);
            [jk jk2]=sort(tmpimg);
            modempr=jk(ceil(numel(jk)*imthr))*1.5;
            %         tmpimg(tmpimg==0)=[];
            %         tmpimg(tmpimg<(median(tmpimg)*.7))=[];
            %         modempr=mode(tmpimg)*3;
            limz.mpr=[0 modempr];
            rgb.mpr=gray(256);
            
            clrs.mprage=rgb.mpr;
            lims.mprage=limz.mpr;
            file.mprage_noskull_afni=[subdir '/mprage.noskull_afni.nii.gz'];
            clrs.mprage_noskull_afni=rgb.mpr;
            lims.mprage_noskull_afni=limz.mpr;
            file.mprage_noskull_fslbet=[subdir '/mprage.noskull_fslbet.nii.gz'];
            clrs.mprage_noskull_fslbet=rgb.mpr;
            lims.mprage_noskull_fslbet=limz.mpr;
            file.mprage_noskull_fs=[subdir '/mprage.noskull_FS.nii.gz'];
            clrs.mprage_noskull_fs=rgb.mpr;
            lims.mprage_noskull_fs=limz.mpr;
            file.mprage_noskull_old=[subdir '/mprage.noskull_old.nii.gz'];
            clrs.mprage_noskull_old=rgb.mpr;
            lims.mprage_noskull_old=limz.mpr;
            file.mprage_noskull=[subdir '/mprage.noskull.nii.gz'];
            clrs.mprage_noskull=rgb.mpr;
            lims.mprage_noskull=limz.mpr;
            file.mprage_noskull_at=[subdir '/mprage.noskull_at.nii.gz'];
            clrs.mprage_noskull_at=rgb.mpr;
            lims.mprage_noskull_at=limz.mpr;
            file.mprage_at=[subdir '/mprage_at.nii.gz'];
            clrs.mprage_at=rgb.mpr;
            lims.mprage_at=limz.mpr;
            file.mprage_noskull_at_EPI=[subdir '/mprage.noskull_at_EPI.nii.gz'];
            clrs.mprage_noskull_at_EPI=rgb.mpr;
            lims.mprage_noskull_at_EPI=limz.mpr;
            file.mprage_at_EPI=[subdir '/mprage_at_EPI.nii.gz'];
            clrs.mprage_at_EPI=rgb.mpr;
            lims.mprage_at_EPI=limz.mpr;
            
            file.mprage_at_afni=[subdir '/mprage_afni_at.nii.gz'];
            clrs.mprage_at=rgb.mpr;
            lims.mprage_at=limz.mpr;
            
            % FS segmentation
            file.aparcaseg=[subdir '/aparc+aseg.nii.gz'];
            clrs.aparcaseg=jet(256);
            lims.aparcaseg=[];
            file.aparcaseg_at=[subdir '/aparc+aseg.atlas.nii.gz'];
            clrs.aparcaseg_at=jet(256);
            lims.aparcaseg_at=[];
            file.aparcaseg_at_EPI=[subdir '/aparc+aseg.atlas_EPI.nii.gz'];
            clrs.aparcaseg_at_EPI=jet(256);
            lims.aparcaseg_at_EPI=[];
            
            % masks of interest
            rgb.WM=[0 0 0;0 128 255]/255; % sky blue for WM masks
            rgb.CSF=[0 0 0;0 0 255]/255; % blue for CSF masks
            rgb.INBRAIN=[0 0 0;255 255 255]/255; % white for in brain
            rgb.GM=[0 0 0;255 0 0]/255; % red for the cortex
            rgb.CBLM=rgb.GM*.8; % darker red for cerebellum
            rgb.SC=rgb.GM*.6; % still darker red for subortex
            
            file.WMero0=[subdir '/aparc+aseg.atlas.WMMASK_ero0.nii.gz'];
            clrs.WMero0=rgb.WM;
            lims.WMero0=[];
            file.WMero1=[subdir '/aparc+aseg.atlas.WMMASK_ero1.nii.gz'];
            clrs.WMero1=rgb.WM;
            lims.WMero1=[];
            file.WMero2=[subdir '/aparc+aseg.atlas.WMMASK_ero2.nii.gz'];
            clrs.WMero2=rgb.WM;
            lims.WMero2=[];
            file.WMero3=[subdir '/aparc+aseg.atlas.WMMASK_ero3.nii.gz'];
            clrs.WMero3=rgb.WM;
            lims.WMero3=[];
            file.WMero4=[subdir '/aparc+aseg.atlas.WMMASK_ero4.nii.gz'];
            clrs.WMero4=rgb.WM;
            lims.WMero4=[];
            file.WMero0_EPI=[subdir '/aparc+aseg.atlas.WMMASK_ero0_EPI.nii.gz'];
            clrs.WMero0_EPI=rgb.WM;
            lims.WMero0_EPI=[];
            file.WMero1_EPI=[subdir '/aparc+aseg.atlas.WMMASK_ero1_EPI.nii.gz'];
            clrs.WMero1_EPI=rgb.WM;
            lims.WMero1_EPI=[];
            file.WMero2_EPI=[subdir '/aparc+aseg.atlas.WMMASK_ero2_EPI.nii.gz'];
            clrs.WMero2_EPI=rgb.WM;
            lims.WMero2_EPI=[];
            file.WMero3_EPI=[subdir '/aparc+aseg.atlas.WMMASK_ero3_EPI.nii.gz'];
            clrs.WMero3_EPI=rgb.WM;
            lims.WMero3_EPI=[];
            file.WMero4_EPI=[subdir '/aparc+aseg.atlas.WMMASK_ero4_EPI.nii.gz'];
            clrs.WMero4_EPI=rgb.WM;
            lims.WMero4_EPI=[];
            file.CSFero0=[subdir '/aparc+aseg.atlas.CSFMASK_ero0.nii.gz'];
            clrs.CSFero0=rgb.CSF;
            lims.CSFero0=[];
            file.CSFero1=[subdir '/aparc+aseg.atlas.CSFMASK_ero1.nii.gz'];
            clrs.CSFero1=rgb.CSF;
            lims.CSFero1=[];
            file.CSFero2=[subdir '/aparc+aseg.atlas.CSFMASK_ero2.nii.gz'];
            clrs.CSFero2=rgb.CSF;
            lims.CSFero2=[];
            file.CSFero3=[subdir '/aparc+aseg.atlas.CSFMASK_ero3.nii.gz'];
            clrs.CSFero3=rgb.CSF;
            lims.CSFero3=[];
            file.CSFero4=[subdir '/aparc+aseg.atlas.CSFMASK_ero4.nii.gz'];
            clrs.CSFero4=rgb.CSF;
            lims.CSFero4=[];
            file.CSFero0_EPI=[subdir '/aparc+aseg.atlas.CSFMASK_ero0_EPI.nii.gz'];
            clrs.CSFero0_EPI=rgb.CSF;
            lims.CSFero0_EPI=[];
            file.CSFero1_EPI=[subdir '/aparc+aseg.atlas.CSFMASK_ero1_EPI.nii.gz'];
            clrs.CSFero1_EPI=rgb.CSF;
            lims.CSFero1_EPI=[];
            file.CSFero2_EPI=[subdir '/aparc+aseg.atlas.CSFMASK_ero2_EPI.nii.gz'];
            clrs.CSFero2_EPI=rgb.CSF;
            lims.CSFero2_EPI=[];
            file.CSFero3_EPI=[subdir '/aparc+aseg.atlas.CSFMASK_ero3_EPI.nii.gz'];
            clrs.CSFero3_EPI=rgb.CSF;
            lims.CSFero3_EPI=[];
            file.CSFero4_EPI=[subdir '/aparc+aseg.atlas.CSFMASK_ero4_EPI.nii.gz'];
            clrs.CSFero4_EPI=rgb.CSF;
            lims.CSFero4_EPI=[];
            file.INBRAINMASK=[subdir '/aparc+aseg.atlas.INBRAINMASK_ero0.nii.gz'];
            clrs.INBRAINMASK=rgb.INBRAIN;
            lims.INBRAINMASK=[];
            file.INBRAINMASK_EPI=[subdir '/aparc+aseg.atlas.INBRAINMASK_ero0_EPI.nii.gz'];
            clrs.INBRAINMASK_EPI=rgb.INBRAIN;
            lims.INBRAINMASK_EPI=[];
            file.GM_RIBBONMASK=[subdir '/aparc+aseg.atlas.GM_RIBBONMASK_ero0.nii.gz'];
            clrs.GM_RIBBONMASK=rgb.GM;
            lims.GM_RIBBONMASK=[];
            file.GM_RIBBONMASK_EPI=[subdir '/aparc+aseg.atlas.GM_RIBBONMASK_ero0_EPI.nii.gz'];
            clrs.GM_RIBBONMASK_EPI=rgb.GM;
            lims.GM_RIBBONMASK_EPI=[];
            file.GM_CBLMMASK_ero0=[subdir '/aparc+aseg.atlas.GM_CBLMMASK_ero0.nii.gz'];
            clrs.GM_CBLMMASK_ero0=rgb.CBLM;
            lims.GM_CBLMMASK_ero0=[];
            file.GM_CBLMMASK_ero0_EPI=[subdir '/aparc+aseg.atlas.GM_CBLMMASK_ero0_EPI.nii.gz'];
            clrs.GM_CBLMMASK_ero0_EPI=rgb.CBLM;
            lims.GM_CBLMMASK_ero0_EPI=[];
            file.GM_CBLMMASK_ero1=[subdir '/aparc+aseg.atlas.GM_CBLMMASK_ero1.nii.gz'];
            clrs.GM_CBLMMASK_ero1=rgb.CBLM;
            lims.GM_CBLMMASK_ero1=[];
            file.GM_CBLMMASK_ero1_EPI=[subdir '/aparc+aseg.atlas.GM_CBLMMASK_ero1_EPI.nii.gz'];
            clrs.GM_CBLMMASK_ero1_EPI=rgb.CBLM;
            lims.GM_CBLMMASK_ero1_EPI=[];
            file.GM_CBLMMASK_ero2=[subdir '/aparc+aseg.atlas.GM_CBLMMASK_ero2.nii.gz'];
            clrs.GM_CBLMMASK_ero2=rgb.CBLM;
            lims.GM_CBLMMASK_ero2=[];
            file.GM_CBLMMASK_ero2_EPI=[subdir '/aparc+aseg.atlas.GM_CBLMMASK_ero2_EPI.nii.gz'];
            clrs.GM_CBLMMASK_ero2_EPI=rgb.CBLM;
            lims.GM_CBLMMASK_ero2_EPI=[];
            file.GM_SCMASK_ero0=[subdir '/aparc+aseg.atlas.GM_SCMASK_ero0.nii.gz'];
            clrs.GM_SCMASK_ero0=rgb.SC;
            lims.GM_SCMASK_ero0=[];
            file.GM_SCMASK_ero0_EPI=[subdir '/aparc+aseg.atlas.GM_SCMASK_ero0_EPI.nii.gz'];
            clrs.GM_SCMASK_ero0_EPI=rgb.SC;
            lims.GM_SCMASK_ero0_EPI=[];
            file.GM_SCMASK_ero1=[subdir '/aparc+aseg.atlas.GM_SCMASK_ero1.nii.gz'];
            clrs.GM_SCMASK_ero1=rgb.SC;
            lims.GM_SCMASK_ero1=[];
            file.GM_SCMASK_ero1_EPI=[subdir '/aparc+aseg.atlas.GM_SCMASK_ero1_EPI.nii.gz'];
            clrs.GM_SCMASK_ero1_EPI=rgb.SC;
            lims.GM_SCMASK_ero1_EPI=[];
            file.GM_SCMASK_ero2=[subdir '/aparc+aseg.atlas.GM_SCMASK_ero2.nii.gz'];
            clrs.GM_SCMASK_ero2=rgb.SC;
            lims.GM_SCMASK_ero2=[];
            file.GM_SCMASK_ero2_EPI=[subdir '/aparc+aseg.atlas.GM_SCMASK_ero2_EPI.nii.gz'];
            clrs.GM_SCMASK_ero2_EPI=rgb.SC;
            lims.GM_SCMASK_ero2_EPI=[];
            file.GM_ALLMASK=[subdir '/aparc+aseg.atlas.GM_ALLMASK_ero0.nii.gz'];
            clrs.GM_ALLMASK=rgb.GM;
            lims.GM_ALLMASK=[];
            file.GM_ALLMASK_EPI=[subdir '/aparc+aseg.atlas.GM_ALLMASK_ero0_EPI.nii.gz'];
            clrs.GM_ALLMASK_EPI=rgb.GM;
            lims.GM_ALLMASK_EPI=[];
            
            iruns=textread([subdir '/runlist.txt']);
            
            for k=iruns
                file.tcat4{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat4.nii.gz'];
                
                % initialize EPI limits
                if (k == iruns(1))
                    tmp=load_untouch_nii(file.tcat4{k});
                    tmpimg=tmp.img(:);
                    [jk jk2]=sort(tmpimg);
                    modeepi=jk(ceil(numel(jk)*imthr))*2;
                    %                 tmpimg(tmpimg==0)=[];
                    %                 tmpimg(tmpimg<(mean(tmpimg)/2))=[];
                    %                 modeepi=mode(tmpimg)*2;
                    limz.tcat=[0 modeepi];
                    rgb.tcat=gray(256);
                end
                file.tcat4_noskull_afni{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat4.noskull_afni.nii.gz'];
                file.tcat4_noskull_fslbet{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat4.noskull_fslbet.nii.gz'];
                file.tcat4_noskull{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat4.noskull.nii.gz'];
                file.tcat4_noskull_at{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat4.noskull_at.nii.gz'];
                file.tcat4_noskull_MPR{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat4.noskull_MPR.nii.gz'];
                file.BOLD{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat.nii.gz'];
                file.BOLD_at{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat.atlas.nii.gz'];
                
                if k>1
                    file.BOLD_at_method2{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat.atlas_motionINT_EPIxEPI1.nii.gz'];
                    file.BOLD_at_method3{k}=[subdir '/pb01.' num2str(sub) '.rest' num2str(iruns(k)) '.tcat.atlas_motionINT_EPIxMPR.nii.gz'];
                end
            end
            clrs.tcat=rgb.tcat;
            lims.tcat=limz.tcat;
            
            % we have now specified all the basic components
            % we can now generate series of images to check registration,
            % segmentation, and other aspects of preprocessing
            ifac=1000000;
            kfac=10;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 0-9 are skull-stripping the mprage
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            imgseq={file.mprage,file.mprage_noskull_afni,file.mprage_noskull_fslbet,file.mprage_noskull_fs,file.mprage_noskull_old,file.mprage_noskull,file.mprage};
            rgbseq={clrs.mprage,clrs.mprage_noskull_afni,clrs.mprage_noskull_fslbet,clrs.mprage_noskull,clrs.mprage_noskull,clrs.mprage_noskull,clrs.mprage};
            limseq={lims.mprage,lims.mprage_noskull_afni,lims.mprage_noskull_fslbet,lims.mprage_noskull_fslbet,lims.mprage_noskull_fslbet,lims.mprage_noskull,lims.mprage};
            olayseq={[1],[2],[3],[4],[5],[6],[7]};
            slices=[.4 .5];
            titleseq={ ...
                {['Subject ' num2str(sub) ],['MPR SKULLSTRIP'],['MP-RAGE']},...
                {['Subject ' num2str(sub) ],['MPR SKULLSTRIP'],['MP-RAGE AFNI skullstrip']},...
                {['Subject ' num2str(sub) ],['MPR SKULLSTRIP'],['MP-RAGE FSL BET skullstrip']},...
                {['Subject ' num2str(sub) ],['MPR SKULLSTRIP'],['MP-RAGE FS skullstrip']},...
                {['Subject ' num2str(sub) ],['MPR SKULLSTRIP'],['MP-RAGE OLD skullstrip']},...
                {['Subject ' num2str(sub) ],['MPR SKULLSTRIP'],['MP-RAGE combined skullstrip']},...
                {['Subject ' num2str(sub) ],['MPR SKULLSTRIP'],['MP-RAGE']},...
                };
            ofile={ ...
                {[picpath '/' num2str(0+ifac*i) '.png']},...
                {[picpath '/' num2str(1+ifac*i) '.png']},...
                {[picpath '/' num2str(2+ifac*i) '.png']},...
                {[picpath '/' num2str(3+ifac*i) '.png']},...
                {[picpath '/' num2str(4+ifac*i) '.png']},...
                {[picpath '/' num2str(5+ifac*i) '.png']},...
                {[picpath '/' num2str(6+ifac*i) '.png']},...
                };
            synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 10-999 are skull-stripping the EPIs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            for k=iruns
                imgseq={file.tcat4{k},file.tcat4_noskull_afni{k},file.tcat4_noskull_fslbet{k},file.tcat4_noskull{k},file.tcat4{k}};
                rgbseq={clrs.tcat,clrs.tcat,clrs.tcat,clrs.tcat,clrs.tcat};
                limseq={lims.tcat,lims.tcat,lims.tcat,lims.tcat,lims.tcat};
                olayseq={[1],[2],[3],[4],[5]};
                slices=[.4 .5 .6];
                titleseq={ ...
                    {['Subject ' num2str(sub) ],['EPI SKULLSTRIP'],['tcat4 rest' num2str(k) ]},...
                    {['Subject ' num2str(sub) ],['EPI SKULLSTRIP'],['tcat4 rest' num2str(k) ' AFNI skullstrip']},...
                    {['Subject ' num2str(sub) ],['EPI SKULLSTRIP'],['tcat4 rest' num2str(k) ' FSL BET skullstrip']},...
                    {['Subject ' num2str(sub) ],['EPI SKULLSTRIP'],['tcat4 rest' num2str(k) ' combined skullstrip']},...
                    {['Subject ' num2str(sub) ],['EPI SKULLSTRIP'],['tcat4 rest' num2str(k) ]},...
                    };
                ofile={ ...
                    {[picpath '/' num2str(0+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(1+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(2+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(3+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(4+(ifac*i)+(kfac*k)) '.png']},...
                    };
                synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 1000-1009 are MPR to ATL registration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offset=1000;
            imgseq={file.TT_N27,file.mprage_noskull_at,file.mprage_at};
            rgbseq={clrs.TT_N27,clrs.mprage_noskull_at,clrs.mprage_at};
            limseq={lims.TT_N27,lims.mprage_noskull_at,lims.mprage_at};
            olayseq={[1],[2],[3]};
            slices=[.4 .5];
            titleseq={ ...
                {['Subject ' num2str(sub) ],['MPR-ATL REGISTRATION'],['MP-RAGE ATLAS (TT-N27)' ]},...
                {['Subject ' num2str(sub) ],['MPR-ATL REGISTRATION'],['MP-RAGE ATLAS-transformed skullstrip']},...
                {['Subject ' num2str(sub) ],['MPR-ATL REGISTRATION'],['MP-RAGE ATLAS-transformed']},...
                };
            ofile={ ...
                {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                };
            synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 1010-1999 are EPI to MPR registration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offset=1000;
            for k=iruns
                imgseq={file.mprage,file.mprage_noskull,file.tcat4_noskull_MPR{k},file.mprage};
                rgbseq={clrs.mprage,clrs.mprage_noskull,clrs.tcat,clrs.mprage};
                limseq={lims.mprage,lims.mprage_noskull,lims.tcat,lims.mprage};
                olayseq={[1],[2],[3],[4]};
                slices=[.4 .5];
                titleseq={ ...
                    {['Subject ' num2str(sub) ],['EPI-MPR REGISTRATION'],['MP-RAGE native space' ]},...
                    {['Subject ' num2str(sub) ],['EPI-MPR REGISTRATION'],['MP-RAGE skullstrip']},...
                    {['Subject ' num2str(sub) ],['EPI-MPR REGISTRATION'],['tcat4 rest' num2str(k) ' skullstripped on MPRAGE']},...
                    {['Subject ' num2str(sub) ],['EPI-MPR REGISTRATION'],['MP-RAGE native space']},...
                    };
                ofile={ ...
                    {[picpath '/' num2str(0+(ifac*i)+offset+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(1+(ifac*i)+offset+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(2+(ifac*i)+offset+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(3+(ifac*i)+offset+(kfac*k)) '.png']},...
                    };
                synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 2000-2999 are EPI to ATL registration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offset=2000;
            
            imgseq={...
                file.mprage_at_EPI,...
                file.GM_RIBBONMASK_EPI,...
                file.GM_CBLMMASK_ero0_EPI,...
                file.GM_SCMASK_ero0_EPI,...
                file.WMero0_EPI...
                file.WMero1_EPI...
                file.WMero2_EPI...
                file.WMero3_EPI...
                file.WMero4_EPI...
                file.CSFero0_EPI...
                file.CSFero1_EPI...
                file.CSFero2_EPI...
                file.CSFero3_EPI...
                file.CSFero4_EPI...
                };
            rgbseq={...
                clrs.mprage_at,...
                clrs.GM_RIBBONMASK,...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                clrs.WMero0...
                clrs.WMero1...
                clrs.WMero2...
                clrs.WMero3...
                clrs.WMero4...
                clrs.CSFero0...
                clrs.CSFero1...
                clrs.CSFero2...
                clrs.CSFero3...
                clrs.CSFero4...
                };
            limseq={...
                lims.mprage,...
                lims.GM_RIBBONMASK,...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                lims.WMero0...
                lims.WMero1...
                lims.WMero2...
                lims.WMero3...
                lims.WMero4...
                lims.CSFero0...
                lims.CSFero1...
                lims.CSFero2...
                lims.CSFero3...
                lims.CSFero4...
                };
            
            
            ctr=numel(limseq)+1;
            imgseq=cat(2,imgseq,file.TT_N27_EPI);
            rgbseq=cat(2,rgbseq,clrs.TT_N27);
            limseq=cat(2,limseq,lims.TT_N27);
            atlctr=ctr;
            ctr=ctr+1;
            imgseq=cat(2,imgseq,file.Power264);
            rgbseq=cat(2,rgbseq,clrs.Power264);
            limseq=cat(2,limseq,lims.Power264);
            roictr=ctr;
            
            
            % start with atlas and mprage
            ofilctr=0;
            olayseq={[atlctr]};
            titleseq={{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas TT-N27 EPI' ]}};
            ofile={{[picpath '/' num2str(0+(ifac*i)+offset) '.png']}};
            ofilctr=ofilctr+1;
            olayseq=cat(2,olayseq,[1]);
            titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas MPRAGE' ]}});
            ofile=cat(2,ofile,{{[picpath '/' num2str(1+(ifac*i)+offset) '.png']}});
            ofilctr=ofilctr+1;
            olayseq=cat(2,olayseq,[1 2 3 4 9 12 ]);
            titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas MPRAGE WM4 CSF2']}});
            ofile=cat(2,ofile,{{[picpath '/' num2str(2+(ifac*i)+offset) '.png']}});
            ofilctr=ofilctr+1;
            olayseq=cat(2,olayseq,[1]);
            titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas MPRAGE' ]}});
            ofile=cat(2,ofile,{{[picpath '/' num2str(3+(ifac*i)+offset) '.png']}});
            
            for k=iruns
                ctr=ctr+1;
                imgseq=cat(2,imgseq,file.BOLD_at{k});
                rgbseq=cat(2,rgbseq,clrs.tcat);
                limseq=cat(2,limseq,lims.tcat);
                
                ofilctr=ofilctr+1;
                olayseq=cat(2,olayseq,[ctr]);
                titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas BOLD Run ' num2str(k) ]}});
                ofile=cat(2,ofile,{{[picpath '/' num2str(0+(ifac*i)+offset+(kfac*k)) '.png']}});
                ofilctr=ofilctr+1;
                olayseq=cat(2,olayseq,[ctr 2 3 4 9 12 ]);
                titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas BOLD Run ' num2str(k) ' WM4 CSF2']}});
                ofile=cat(2,ofile,{{[picpath '/' num2str(1+(ifac*i)+offset+(kfac*k)) '.png']}});
                ofilctr=ofilctr+1;
                olayseq=cat(2,olayseq,[ctr roictr ]);
                titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas BOLD Run ' num2str(k) ' Power264']}});
                ofile=cat(2,ofile,{{[picpath '/' num2str(2+(ifac*i)+offset+(kfac*k)) '.png']}});
                ofilctr=ofilctr+1;
                olayseq=cat(2,olayseq,[ctr]);
                titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas BOLD Run ' num2str(k) ]}});
                ofile=cat(2,ofile,{{[picpath '/' num2str(3+(ifac*i)+offset+(kfac*k)) '.png']}});
            end
            
            slices=[.4 .5 .6];
            synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 3000-3999 are mask images on MPRAGEs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offset=3000;
            imgseq={...
                file.mprage_at,...
                file.GM_RIBBONMASK,...
                file.WMero0...
                file.WMero1...
                file.WMero2...
                file.WMero3...
                file.WMero4...
                file.CSFero0...
                file.CSFero1...
                file.CSFero2...
                file.CSFero3...
                file.CSFero4...
                file.GM_CBLMMASK_ero0,...
                file.GM_CBLMMASK_ero1,...
                file.GM_CBLMMASK_ero2,...
                file.GM_SCMASK_ero0,...
                file.GM_SCMASK_ero1,...
                file.GM_SCMASK_ero2,...
                };
            rgbseq={...
                clrs.mprage_at,...
                clrs.GM_RIBBONMASK,...
                clrs.WMero0...
                clrs.WMero1...
                clrs.WMero2...
                clrs.WMero3...
                clrs.WMero4...
                clrs.CSFero0...
                clrs.CSFero1...
                clrs.CSFero2...
                clrs.CSFero3...
                clrs.CSFero4...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                };
            limseq={...
                lims.mprage,...
                lims.GM_RIBBONMASK,...
                lims.WMero0...
                lims.WMero1...
                lims.WMero2...
                lims.WMero3...
                lims.WMero4...
                lims.CSFero0...
                lims.CSFero1...
                lims.CSFero2...
                lims.CSFero3...
                lims.CSFero4...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                };
            olayseq={[1],[1 2 13 16 3 8],[1 2 13 16 4 9],[1 2 13 16 5 10],[1 2 13 16 6 11],[1 2 13 16 7 12]};
            slices=[.4 .5];
            titleseq={ ...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE' ]},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero0']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero1']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero2']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero3']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero4']},...
                };
            ofile={ ...
                {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(3+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(4+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(5+(ifac*i)+offset) '.png']},...
                };
            synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 4000-4999 are mask images on MPRAGEs at EPI resolution
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offset=4000;
            imgseq={...
                file.mprage_at_EPI,...
                file.GM_RIBBONMASK_EPI,...
                file.GM_CBLMMASK_ero0_EPI,...
                file.GM_SCMASK_ero0_EPI,...
                file.WMero0_EPI...
                file.WMero1_EPI...
                file.WMero2_EPI...
                file.WMero3_EPI...
                file.WMero4_EPI...
                file.CSFero0_EPI...
                file.CSFero1_EPI...
                file.CSFero2_EPI...
                file.CSFero3_EPI...
                file.CSFero4_EPI...
                file.GM_CBLMMASK_ero0_EPI,...
                file.GM_CBLMMASK_ero1_EPI,...
                file.GM_CBLMMASK_ero2_EPI,...
                file.GM_SCMASK_ero0_EPI,...
                file.GM_SCMASK_ero1_EPI,...
                file.GM_SCMASK_ero2_EPI,...
                };
            rgbseq={...
                clrs.mprage_at,...
                clrs.GM_RIBBONMASK,...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                clrs.WMero0...
                clrs.WMero1...
                clrs.WMero2...
                clrs.WMero3...
                clrs.WMero4...
                clrs.CSFero0...
                clrs.CSFero1...
                clrs.CSFero2...
                clrs.CSFero3...
                clrs.CSFero4...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                };
            limseq={...
                lims.mprage,...
                lims.GM_RIBBONMASK,...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                lims.WMero0...
                lims.WMero1...
                lims.WMero2...
                lims.WMero3...
                lims.WMero4...
                lims.CSFero0...
                lims.CSFero1...
                lims.CSFero2...
                lims.CSFero3...
                lims.CSFero4...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                };
            olayseq={[1],[1 2 15 18 5 10],[1 2 15 18 6 11],[1 2 15 18 7 12],[1 2 15 18 8 13],[1 2 15 18 9 14],[1],[1 2 3 4 9 12]};
            slices=[.4 .5 .6];
            titleseq={ ...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE' ]},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero0']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero1']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero2']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero3']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments ero4']},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE' ]},...
                {['Subject ' num2str(sub) ],['NUISANCE MASK EROSION'],['MP-RAGE + FS segments WM4 CSF2']},...
                };
            ofile={ ...
                {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(3+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(4+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(5+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(6+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(7+(ifac*i)+offset) '.png']},...
                };
            synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 5000-5999 are verions of EPI to ATL registration
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offset=5000;
            
            imgseq={...
                file.mprage_at_EPI,...
                file.GM_RIBBONMASK_EPI,...
                file.GM_CBLMMASK_ero0_EPI,...
                file.GM_SCMASK_ero0_EPI,...
                file.WMero0_EPI...
                file.WMero1_EPI...
                file.WMero2_EPI...
                file.WMero3_EPI...
                file.WMero4_EPI...
                file.CSFero0_EPI...
                file.CSFero1_EPI...
                file.CSFero2_EPI...
                file.CSFero3_EPI...
                file.CSFero4_EPI...
                };
            rgbseq={...
                clrs.mprage_at,...
                clrs.GM_RIBBONMASK,...
                clrs.GM_CBLMMASK_ero0,...
                clrs.GM_SCMASK_ero0,...
                clrs.WMero0...
                clrs.WMero1...
                clrs.WMero2...
                clrs.WMero3...
                clrs.WMero4...
                clrs.CSFero0...
                clrs.CSFero1...
                clrs.CSFero2...
                clrs.CSFero3...
                clrs.CSFero4...
                };
            limseq={...
                lims.mprage,...
                lims.GM_RIBBONMASK,...
                lims.GM_CBLMMASK_ero0,...
                lims.GM_SCMASK_ero0,...
                lims.WMero0...
                lims.WMero1...
                lims.WMero2...
                lims.WMero3...
                lims.WMero4...
                lims.CSFero0...
                lims.CSFero1...
                lims.CSFero2...
                lims.CSFero3...
                lims.CSFero4...
                };
            
            
            ctr=numel(limseq)+1;
            imgseq=cat(2,imgseq,file.TT_N27_EPI);
            rgbseq=cat(2,rgbseq,clrs.TT_N27);
            limseq=cat(2,limseq,lims.TT_N27);
            atlctr=ctr;
            ctr=ctr+1;
            imgseq=cat(2,imgseq,file.Power264);
            rgbseq=cat(2,rgbseq,clrs.Power264);
            limseq=cat(2,limseq,lims.Power264);
            roictr=ctr;
            
            
            % start with atlas and mprage
            ofilctr=0;
            olayseq={[atlctr]};
            titleseq={{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION'],['Atlas TT-N27 EPI' ]}};
            ofile={{[picpath '/' num2str(0+(ifac*i)+offset) '.png']}};
            ofilctr=ofilctr+1;
            olayseq=cat(2,olayseq,[1]);
            titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION METHOD1'],['Atlas MPRAGE' ]}});
            ofile=cat(2,ofile,{{[picpath '/' num2str(1+(ifac*i)+offset) '.png']}});
            ofilctr=ofilctr+1;
            olayseq=cat(2,olayseq,[1 2 3 4 9 12 ]);
            titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION METHOD1'],['Atlas MPRAGE WM4 CSF2']}});
            ofile=cat(2,ofile,{{[picpath '/' num2str(2+(ifac*i)+offset) '.png']}});
            ofilctr=ofilctr+1;
            olayseq=cat(2,olayseq,[1]);
            titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION METHOD1'],['Atlas MPRAGE' ]}});
            ofile=cat(2,ofile,{{[picpath '/' num2str(3+(ifac*i)+offset) '.png']}});
            
            %         runctr=0;
            %         for k=iruns
            %             runctr=runctr+1;
            %             ctr=ctr+1;
            %             imgseq=cat(2,imgseq,file.BOLD{k});
            %             rgbseq=cat(2,rgbseq,clrs.tcat);
            %             limseq=cat(2,limseq,lims.tcat);
            %             ftypeseq=[ftypeseq 5];
            %
            %             ofilctr=ofilctr+1;
            %             olayseq=cat(2,olayseq,[ctr]);
            %             titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION NATIVE'],['BOLD Run ' num2str(k) ]}});
            %             ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+10) '.png']}});
            %         end
            
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                ctr=ctr+1;
                imgseq=cat(2,imgseq,file.BOLD_at{k});
                rgbseq=cat(2,rgbseq,clrs.tcat);
                limseq=cat(2,limseq,lims.tcat);
                
                ofilctr=ofilctr+1;
                olayseq=cat(2,olayseq,[ctr]);
                titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION METHOD1'],['MOT-EPI1,EPI1-MPR,MPR-ATL'],['Atlas BOLD Run ' num2str(k) ]}});
                ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+100) '.png']}});
            end
            
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                ctr=ctr+1;
                if runctr==1
                    imgseq=cat(2,imgseq,file.BOLD_at{k});
                else
                    imgseq=cat(2,imgseq,file.BOLD_at_method2{k});
                end
                rgbseq=cat(2,rgbseq,clrs.tcat);
                limseq=cat(2,limseq,lims.tcat);
                
                ofilctr=ofilctr+1;
                olayseq=cat(2,olayseq,[ctr]);
                titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION METHOD2'],['MOT-EPIX,EPIX-EPI1,EPI1-MPR,MPR-ATL'],['Atlas BOLD Run ' num2str(k) ]}});
                ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+200) '.png']}});
            end
            
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                ctr=ctr+1;
                if runctr==1
                    imgseq=cat(2,imgseq,file.BOLD_at{k});
                else
                    imgseq=cat(2,imgseq,file.BOLD_at_method3{k});
                end
                rgbseq=cat(2,rgbseq,clrs.tcat);
                limseq=cat(2,limseq,lims.tcat);
                
                ofilctr=ofilctr+1;
                olayseq=cat(2,olayseq,[ctr]);
                titleseq=cat(2,titleseq,{{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION METHOD3'],['MOT-EPIX,EPIX-MPR,MPR-ATL'],['Atlas BOLD Run ' num2str(k) ]}});
                ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+300) '.png']}});
            end
            
            slices=[.4 .5 .6];
            synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            
            % different resolution probably, have to do them separately
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                imgseq={file.BOLD{k}};
                rgbseq={clrs.tcat};
                limseq={lims.tcat};
                
                olayseq={[1]};
                titleseq={{['Subject ' num2str(sub) ],['EPI-ATL REGISTRATION NATIVE'],['BOLD Run ' num2str(k) ]}};
                ofile={{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+10) '.png']}};
                
                slices=[.4 .5 .6];
                synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            end
            
            
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% 6000-6009 are MPR-ATL under different skull strips
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            offset=6000;
            for k=iruns
                imgseq={file.TT_N27,file.mprage_noskull_at,file.mprage_at,file.mprage_at_afni};
                rgbseq={clrs.TT_N27,clrs.mprage_noskull_at,clrs.mprage_at,clrs.mprage_at};
                limseq={lims.TT_N27,lims.mprage_noskull_at,lims.mprage_at,lims.mprage_at};
                olayseq={[1],[2],[3],[4],[1]};
                slices=[.4 .5];
                titleseq={ ...
                    {['Subject ' num2str(sub) ],['EFFECT OF SKULLSTRIP'],['MP-RAGE ATLAS (TT-N27)' ]},...
                    {['Subject ' num2str(sub) ],['EFFECT OF SKULLSTRIP'],['MP-RAGE ATLAS-transformed skullstrip']},...
                    {['Subject ' num2str(sub) ],['EFFECT OF SKULLSTRIP'],['MP-RAGE ATLAS-transformed']},...
                    {['Subject ' num2str(sub) ],['EFFECT OF SKULLSTRIP'],['MP-RAGE ATLAS-transformed using AFNI skullstrip']},...
                    {['Subject ' num2str(sub) ],['EFFECT OF SKULLSTRIP'],['MP-RAGE ATLAS (TT-N27)' ]},...
                    };
                ofile={ ...
                    {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                    {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                    {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                    {[picpath '/' num2str(3+(ifac*i)+offset) '.png']},...
                    {[picpath '/' num2str(4+(ifac*i)+offset) '.png']},...
                    };
                synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile);
            end
            
        end
        
        
    case 'assemblepictures'
        
        % turns out it was too much information the first time around, group
        % the pictures into more manageable chunks
        
        subsdir=[pdir '/subs'];
        basepicpath=[logpath '/SUMMARY_pics'];
        %         picpath=[basepicpath '/rawpics'];
        %         system(['mkdir -p ' picpath ]);
        
        switch scriptstr
            case {'JDP_preprocess.sh','JDP_P1_preprocess.sh'}
                resdir='P1.results';
        end
        
        ifac=1000000;
        kfac=10;
        
        % folder for skullstripping the mprages
        fprintf('Copying MPRAGE skullstrip images\n');
        subfolder=[basepicpath '/MPRnoskull'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            ofile={ ...
                {[picpath '/' num2str(0+ifac*i) '.png']},...
                {[picpath '/' num2str(1+ifac*i) '.png']},...
                {[picpath '/' num2str(2+ifac*i) '.png']},...
                {[picpath '/' num2str(3+ifac*i) '.png']},...
                {[picpath '/' num2str(4+ifac*i) '.png']},...
                {[picpath '/' num2str(5+ifac*i) '.png']},...
                {[picpath '/' num2str(6+ifac*i) '.png']},...
                };
            for j=1:numel(ofile)
                str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
            end
        end
        
        % folder for mprage to atlas transformation
        fprintf('Copying MPRAGE atlas images\n');
        offset=1000;
        subfolder=[basepicpath '/MPRatlas'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            ofile={ ...
                {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                };
            for j=1:numel(ofile)
                str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
            end
        end
        
        % folder for segmentations on the mprage
        fprintf('Copying MPRAGE segmentation images\n');
        offset=3000;
        subfolder=[basepicpath '/MPRsegment'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            ofile={ ...
                {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(3+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(4+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(5+(ifac*i)+offset) '.png']},...
                };
            for j=1:numel(ofile)
                str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
            end
        end
        
        % folder for segmentations on the mprage at EPI resolution
        fprintf('Copying MPRAGE EPI images\n');
        offset=4000;
        subfolder=[basepicpath '/MPRsegmentEPI'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            ofile={ ...
                {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(3+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(4+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(5+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(6+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(7+(ifac*i)+offset) '.png']},...
                };
            for j=1:numel(ofile)
                str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
            end
        end
        
        % folder for skullstripping the EPI
        fprintf('Copying EPI skullstrip images\n');
        subfolder=[basepicpath '/EPInoskull'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            iruns=textread([subdir '/runlist.txt']);
            for k=iruns
                ofile={ ...
                    {[picpath '/' num2str(0+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(1+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(2+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(3+(ifac*i)+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(4+(ifac*i)+(kfac*k)) '.png']},...
                    };
                for j=1:numel(ofile)
                    str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
                end
                
            end
        end
        
        % folder for registering EPI to the MPRAGE
        fprintf('Copying EPI MPRAGE images\n');
        offset=1000;
        subfolder=[basepicpath '/EPIonMPR'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            iruns=textread([subdir '/runlist.txt']);
            for k=iruns
                ofile={ ...
                    {[picpath '/' num2str(1+(ifac*i)+offset+(kfac*k)) '.png']},...
                    {[picpath '/' num2str(2+(ifac*i)+offset+(kfac*k)) '.png']},...
                    };
                for j=1:numel(ofile)
                    str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
                end
                
            end
        end
        
        % folder for EPI segmentation in atlas space
        fprintf('Copying EPI atlas images\n');
        subfolder=[basepicpath '/EPIatlas_segment'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            iruns=textread([subdir '/runlist.txt']);
            
            offset=2000;
            
            % start with atlas and mprage
            ofilctr=0;
            ofile={{[picpath '/' num2str(0+(ifac*i)+offset) '.png']}};
            ofilctr=ofilctr+1;
            ofile=cat(2,ofile,{{[picpath '/' num2str(1+(ifac*i)+offset) '.png']}});
            ofilctr=ofilctr+1;
            ofile=cat(2,ofile,{{[picpath '/' num2str(2+(ifac*i)+offset) '.png']}});
            ofilctr=ofilctr+1;
            ofile=cat(2,ofile,{{[picpath '/' num2str(3+(ifac*i)+offset) '.png']}});
            
            for k=iruns
                
                ofilctr=ofilctr+1;
                ofile=cat(2,ofile,{{[picpath '/' num2str(0+(ifac*i)+offset+(kfac*k)) '.png']}});
                ofilctr=ofilctr+1;
                ofile=cat(2,ofile,{{[picpath '/' num2str(1+(ifac*i)+offset+(kfac*k)) '.png']}});
                ofilctr=ofilctr+1;
                ofile=cat(2,ofile,{{[picpath '/' num2str(2+(ifac*i)+offset+(kfac*k)) '.png']}});
                ofilctr=ofilctr+1;
                ofile=cat(2,ofile,{{[picpath '/' num2str(3+(ifac*i)+offset+(kfac*k)) '.png']}});
            end
            
            for j=1:numel(ofile)
                str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
            end
            
            
        end
        
        % folder for MPR-ATL via AFNI- vs consensus-skullstrip
        fprintf('Copying methods for registering MPRAGEs\n');
        offset=6000;
        subfolder=[basepicpath '/MPRatlas_methods'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            ofile={ ...
                {[picpath '/' num2str(0+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(1+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(2+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(3+(ifac*i)+offset) '.png']},...
                {[picpath '/' num2str(4+(ifac*i)+offset) '.png']},...
                };
            for j=1:numel(ofile)
                str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
            end
        end
        clear ofile;
        
        
        % folder for various EPI-ATL registration methods
        fprintf('Copying methods for EPI atlas registration\n');
        subfolder=[basepicpath '/EPIatlas_methods'];
        str=['rm -rf ' subfolder]; system(str);
        str=['mkdir -p ' subfolder]; system(str);
        for i=startsub:endsub
            
            sub=subjlist{i};
            subdir=[subsdir '/' num2str(sub) '/' resdir];
            picpath=[subsdir '/' num2str(sub) '/P1.pictures'];
            iruns=textread([subdir '/runlist.txt']);
            
            offset=5000;
            
            % start with atlas and mprage
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                if k==iruns(1)
                    ofile={{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+10) '.png']}};
                else
                    ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+10) '.png']}});
                end
            end
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+100) '.png']}});
            end
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+200) '.png']}});
            end
            runctr=0;
            for k=iruns
                runctr=runctr+1;
                ofile=cat(2,ofile,{{[picpath '/' num2str((runctr-1)+(ifac*i)+offset+300) '.png']}});
            end
            
            for j=1:numel(ofile)
                str=['cp ' ofile{j}{1} ' ' subfolder ]; system(str);
            end
            
            
        end
        
        ofile=[logpath '/JDP_P1_preprocess_DATASELECT.txt'];
        fid=fopen(ofile,'w');
        for i=startsub:endsub
            sub=subjlist{i};
            subdir=[datadir '/' num2str(sub) '/' resdir];
            iruns=textread([subdir '/runlist.txt']);
            for k=iruns
                fprintf(fid,'%s\t%d\t1\n',sub,iruns(k));
            end
            
        end
        fclose(fid);
        
        
    case 'makemats'

        
        % for each subject, gather relevant information and write it to a .mat
        for i=startsub:endsub
            fprintf('Creating .mat file for subject %d of %d through %d\n',i,startsub,endsub);
            clear QC;
            QC=giveQCmat(scriptdir,datadir,subjlist,i);
            QCname=[datadir '/' QC.sub '/P1.mat' ];
            system(['rm -rf ' QCname ]);
            save(QCname,'QC','-v7.3');
        end
        
      
        
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [img]=showbrik(V,orien,slc,t,varargin)

%d1ATL: left to right; front down; bottom left
%d2ATL: back to front; left top; bottom left
%d3ATL: bottom to top; left top; front right
%d1EPI: left to right; front down; bottom left
%d2EPI: back to front; left top; bottom left
%d3EPI: bottom to top; left top; front right

% presumes RAS orientation, just saying.

switch orien
    case 'side' % rot 90 ccw
        img=(rot90(squeeze(V(slc,:,:,t))));
    case 'front' % rot 90 ccw
        img=(rot90(squeeze(V(:,slc,:,t))));
    case 'top' % rot 90 ccw
        img=(rot90(squeeze(V(:,:,slc,t))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%I
function [ming]=multslc(V,orien,slcz,t,varargin)

ming=[];

for z=1:numel(slcz)
    [img]=showbrik(V,orien,slcz(z),t);
    if isempty(varargin)
        ming=[ming img];
    else
        ming=[ming; img];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rgbval] = convert2rgb(tmp,scl,varargin)

% tmp is an image
% scl is like jet or copper or self-made
% varargin is limits for scaling the image



%if varargin then we're scaling image
if ~isempty(varargin)
    minz=varargin{1,1}(1);
    maxz=varargin{1,1}(2);
    
    % add an extra column with limits
    tmp=[tmp zeros(size(tmp,1),1)];
    tmp(end,end)=minz;
    tmp(end-1,end)=maxz;
    
    tmp(tmp<minz)=minz;
    tmp(tmp>maxz)=maxz;
    tmp=tmp-min(tmp(:));
    tmp=tmp/max(tmp(:));
    tmp=uint8(tmp*256);
    
    % delete the extra column
    tmp(:,end)=[];
else
    tmp=tmp-min(tmp(:));
    tmp=tmp/max(tmp(:));
    tmp=uint8(tmp*256);
end
rgbval=ind2rgb(tmp,scl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nimg] = olay(bimg,timg,msk)

% overlay one image on another

msk=~~msk;
msk=repmat(msk,[1 1 size(bimg,3)]);
nimg=bimg;
nimg(msk)=timg(msk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function synthesize_img(imgseq,rgbseq,limseq,olayseq,slices,titleseq,ofile)

% load in the images, they will have similar orientations and sizes
d=numel(imgseq);
for i=1:d
   
    % open the nii, convert to RAS, return pixel dimensions
    % note that pixel dims need to be the same for all images
    [data{i} dims]=nii_give_RAS(imgseq{i});
end

% these are the slices we'll illustrate
dd=size(data{1});
sliceind=ceil(dd'*slices);

% extract the slices in 3 orientations, convert to rgb values
for i=1:d
    img{i,1}=double(multslc(data{i},'side',sliceind(1,:),1,1));
    img{i,2}=double(multslc(data{i},'front',sliceind(2,:),1,1));
    img{i,3}=double(multslc(data{i},'top',sliceind(3,:),1,1));
    if isempty(limseq{i})
        limseq{i}=[min(data{i}(:)) max(data{i}(:))];
    end
    rgb{i,1}=convert2rgb(img{i,1},rgbseq{i},limseq{i});
    rgb{i,2}=convert2rgb(img{i,2},rgbseq{i},limseq{i});
    rgb{i,3}=convert2rgb(img{i,3},rgbseq{i},limseq{i});
    
end

% perform any desired overlays of the images
for i=1:numel(olayseq)
    if numel(olayseq{i})==1
        im{i,1}=rgb{olayseq{i},1};
        im{i,2}=rgb{olayseq{i},2};
        im{i,3}=rgb{olayseq{i},3};
    else
        for k=1:numel(olayseq{i})
            if k==1
                im{i,1}=rgb{olayseq{i}(k),1};
                im{i,2}=rgb{olayseq{i}(k),2};
                im{i,3}=rgb{olayseq{i}(k),3};
            else
                im{i,1}=olay(im{i,1},rgb{olayseq{i}(k),1},img{olayseq{i}(k),1});
                im{i,2}=olay(im{i,2},rgb{olayseq{i}(k),2},img{olayseq{i}(k),2});
                im{i,3}=olay(im{i,3},rgb{olayseq{i}(k),3},img{olayseq{i}(k),3});
            end
        end
    end
end

for i=1:numel(olayseq)
    % write out the pictures
    close all;
    h=figure;
    set(h,'position',[10 10 1920 1080]);
    set(h,'visible','off');
    fs=16;
    
    % 1, 2, 3 are side, front, and top
    % if anisotropic need to scale with daspect
    
    subplot(1,3,1);
    image(im{i,1});
    daspect([dims(3) dims(2) dims(1)]); %z,y
    axis('off');
    subplot(1,3,2);
    image(im{i,2});
    axis('off');
    daspect([dims(3) dims(1) dims(2)]); %z,z
    hh=title(titleseq{i},'Interpreter','none'); set(hh,'fontsize',fs);
    subplot(1,3,3);
    image(im{i,3});
    daspect([dims(2) dims(1) dims(3)]); %y,x
    axis('off');
    
    % for some reason the white background looks bad to me
    whitebg(h,[.01 .01 .01]);
    set(gcf,'InvertHardcopy','off');
    
    
    set(h,'paperpositionmode','auto');
    str=ofile{i};
    print(gcf,'-dpng',str{1});
    
    close(h);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QC2] = giveQCmat(scriptdir,datadir,isub,i)


QC(i).sub=isub{i};
QC(i).subdatadir=[datadir '/' QC(i).sub '/P1.results' ];

% gather the structural images
tmp=nii_give_RAS([QC(i).subdatadir '/mprage.nii.gz']);
QC(i).mprage.orig=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/mprage.noskull.nii.gz']);
QC(i).mprage.noskull=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/mprage_at.nii.gz']);
QC(i).mprage.at=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/mprage_at_EPI.nii.gz']);
QC(i).mprage.at_EPI=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/mprage.noskull_at.nii.gz']);
QC(i).mprage.noskull_at=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/mprage.noskull_at_EPI.nii.gz']);
QC(i).mprage.noskull_at_EPI=single(tmp);

% gather atlas images
tmp=nii_give_RAS([QC(i).subdatadir '/TT_N27.nii.gz']);
QC(i).atlas.TT_N27=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/TT_N27_EPI.nii.gz']);
QC(i).atlas.TT_N27_EPI=single(tmp);

% gather the mask images
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.nii.gz']);
QC(i).aparcaseg.orig=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.nii.gz']);
QC(i).aparcaseg.at=single(tmp);
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas_EPI.nii.gz']);
QC(i).aparcaseg.at_EPI=single(tmp);

tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero0.nii.gz']);
QC(i).WMMASK.ero0=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero1.nii.gz']);
QC(i).WMMASK.ero1=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero2.nii.gz']);
QC(i).WMMASK.ero2=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero3.nii.gz']);
QC(i).WMMASK.ero3=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero4.nii.gz']);
QC(i).WMMASK.ero4=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero0_EPI.nii.gz']);
QC(i).WMMASK.ero0_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero1_EPI.nii.gz']);
QC(i).WMMASK.ero1_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero2_EPI.nii.gz']);
QC(i).WMMASK.ero2_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero3_EPI.nii.gz']);
QC(i).WMMASK.ero3_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.WMMASK_ero4_EPI.nii.gz']);
QC(i).WMMASK.ero4_EPI=~~tmp;

tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero0.nii.gz']);
QC(i).CSFMASK.ero0=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero1.nii.gz']);
QC(i).CSFMASK.ero1=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero2.nii.gz']);
QC(i).CSFMASK.ero2=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero3.nii.gz']);
QC(i).CSFMASK.ero3=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero4.nii.gz']);
QC(i).CSFMASK.ero4=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero0_EPI.nii.gz']);
QC(i).CSFMASK.ero0_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero1_EPI.nii.gz']);
QC(i).CSFMASK.ero1_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero2_EPI.nii.gz']);
QC(i).CSFMASK.ero2_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero3_EPI.nii.gz']);
QC(i).CSFMASK.ero3_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.CSFMASK_ero4_EPI.nii.gz']);
QC(i).CSFMASK.ero4_EPI=~~tmp;

tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_ALLMASK_ero0.nii.gz']);
QC(i).GMMASK.ALLMASK_ero0=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_ALLMASK_ero0_EPI.nii.gz']);
QC(i).GMMASK.ALLMASK_ero0_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_RIBBONMASK_ero0.nii.gz']);
QC(i).GMMASK.RIBBONMASK_ero0=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_RIBBONMASK_ero0_EPI.nii.gz']);
QC(i).GMMASK.RIBBONMASK_ero0_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.INBRAINMASK_ero0.nii.gz']);
QC(i).INBRAINMASK.ero0=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.INBRAINMASK_ero0_EPI.nii.gz']);
QC(i).INBRAINMASK.ero0_EPI=~~tmp;

tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_SCMASK_ero0.nii.gz']);
QC(i).GMMASK.SCMASK_ero0=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_SCMASK_ero1.nii.gz']);
QC(i).GMMASK.SCMASK_ero1=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_SCMASK_ero2.nii.gz']);
QC(i).GMMASK.SCMASK_ero2=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_SCMASK_ero0_EPI.nii.gz']);
QC(i).GMMASK.SCMASK_ero0_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_SCMASK_ero1_EPI.nii.gz']);
QC(i).GMMASK.SCMASK_ero1_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_SCMASK_ero2_EPI.nii.gz']);
QC(i).GMMASK.SCMASK_ero2_EPI=~~tmp;

tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_CBLMMASK_ero0.nii.gz']);
QC(i).GMMASK.CBLMMASK_ero0=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_CBLMMASK_ero1.nii.gz']);
QC(i).GMMASK.CBLMMASK_ero1=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_CBLMMASK_ero2.nii.gz']);
QC(i).GMMASK.CBLMMASK_ero2=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_CBLMMASK_ero0_EPI.nii.gz']);
QC(i).GMMASK.CBLMMASK_ero0_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_CBLMMASK_ero1_EPI.nii.gz']);
QC(i).GMMASK.CBLMMASK_ero1_EPI=~~tmp;
tmp=nii_give_RAS([QC(i).subdatadir '/aparc+aseg.atlas.GM_CBLMMASK_ero2_EPI.nii.gz']);
QC(i).GMMASK.CBLMMASK_ero2_EPI=~~tmp;

% gather run information
QC(i).iruns=textread([QC(i).subdatadir '/runlist.txt']);
[QC(i).iruns_convert{1} QC(i).iruns_convert{2}]=textread([QC(i).subdatadir '/runconvertfile.txt'],'%s%s');

for j=QC(i).iruns
    
    % grab the TR and slicetime and skiptrs inforamtion
    for k=1:numel(QC(i).iruns_convert{2})
        [jk nm jk1]=fileparts(QC(i).iruns_convert{2}{k});
        [jk nm jk1]=fileparts(nm);
        QC(i).TR{k}=textread([QC(i).subdatadir '/' nm '.TR']);
        QC(i).SLICETIME{k}=textread([QC(i).subdatadir '/' nm '.SLICETIME'],'%s');
        QC(i).SKIPTRS{k}=textread([QC(i).subdatadir '/' nm '.SKIPTRS'],'%s');
    end
    
    
    
    % extract features of the EPI data
    %
    % this command will give the relevant images
    % str = ['ls ' QC(i).subdatadir '/*.nii.gz | grep -v aparc | grep -v mprage | grep atlas | grep rest' num2str(j)];
    QC(i).imglist{j}= {...
        [QC(i).subdatadir '/pb01.' QC(i).sub '.rest' num2str(j) '.tcat.atlas.nii.gz'],...
        [QC(i).subdatadir '/pb02.' QC(i).sub '.rest' num2str(j) '.despike.atlas.nii.gz'],...
        [QC(i).subdatadir '/pb03.' QC(i).sub '.rest' num2str(j) '.despike.tshift.atlas.nii.gz'],...
        [QC(i).subdatadir '/pb04.' QC(i).sub '.rest' num2str(j) '.tshift.atlas.nii.gz'],...
        [QC(i).subdatadir '/pb05.' QC(i).sub '.rest' num2str(j) '.tshift.despike.atlas.nii.gz'],...
        [QC(i).subdatadir '/pb06.' QC(i).sub '.rest' num2str(j) '.tcat.atlas.tshift.nii.gz'],...
        [QC(i).subdatadir '/pb07.' QC(i).sub '.rest' num2str(j) '.despike.atlas.tshift.nii.gz'],...
        };
    
    % extract pertinent features of each run
    for k=1:numel(QC(i).imglist{j})
        
        % load the image
        tmp=nii_give_RAS(QC(i).imglist{j}{k});
        
        % calculate the mean (ANATAVE)
        QC(i).EPI.ANATAVE{j,k}=nanmean(tmp,4);
        QC(i).EPI.ANATAVE_STD{j,k}=nanstd(double(tmp),[],4);
        
        % extract representative GM, WM, and CSF voxels
        % this is just for quick visuals:
        % - avoid susceptibility zones (ANATAVE>200)
        % - trim to a smallish number of randomly selected voxels
        % - within some FS-defined tissue compartment
        d=size(tmp);
        flattmp=double(reshape(tmp,[d(1)*d(2)*d(3) d(4)]));
        
        % choose the voxels on the first image, keep the same after
        if k==1
            tmpthr=max(QC(i).EPI.ANATAVE{j,k}(:))*0.2;
            susMASK=QC(i).EPI.ANATAVE{j,k}(:)>tmpthr;
            
            tmpMASK=QC(i).CSFMASK.ero2_EPI(:);
            tmpMASK=susMASK & tmpMASK;
            remnum=100;
            if nnz(tmpMASK)>remnum
                findrems=find(tmpMASK);
                [jk in]=sort(rand(size(findrems,1),1));
                findrems=findrems(in);
                settozero=findrems(remnum+1:end);
                tmpMASK(settozero)=0;
            end
            tmpCSFMASK=tmpMASK;
            
            tmpMASK=QC(i).WMMASK.ero4_EPI(:);
            tmpMASK=susMASK & tmpMASK;
            remnum=500;
            if nnz(tmpMASK)>remnum
                findrems=find(tmpMASK);
                [jk in]=sort(rand(size(findrems,1),1));
                findrems=findrems(in);
                settozero=findrems(remnum+1:end);
                tmpMASK(settozero)=0;
            end
            tmpWMMASK=tmpMASK;
            
            tmpMASK=QC(i).GMMASK.ALLMASK_ero0_EPI(:);
            tmpMASK=susMASK & tmpMASK;
            remnum=1000;
            if nnz(tmpMASK)>remnum
                findrems=find(tmpMASK);
                [jk in]=sort(rand(size(findrems,1),1));
                findrems=findrems(in);
                settozero=findrems(remnum+1:end);
                tmpMASK(settozero)=0;
            end
            tmpGMMASK=tmpMASK;
        end
        
        QC(i).EPI.WM{j,k}=single(flattmp(tmpWMMASK,:));
        QC(i).EPI.CSF{j,k}=single(flattmp(tmpCSFMASK,:));
        QC(i).EPI.GM{j,k}=single(flattmp(tmpGMMASK,:));
        
        % extract the 264 timeseries
        power264=[scriptdir '/Power264_indiv.nii.gz'];
        QC(i).EPI.Power264{j,k}=give264tcs(flattmp,power264);
        
        
        % obtain MEAN, SD, and DVARS calculations
        QC(i).EPI.MEAN{j,k}=nanmean(flattmp(QC(i).INBRAINMASK.ero0_EPI(:),:),1);
        QC(i).EPI.SD{j,k}=nanstd(double(flattmp(QC(i).INBRAINMASK.ero0_EPI(:),:)),[],1);
        tmp2=flattmp;
        tmp2=flattmp-repmat(nanmean(flattmp,2),[1 size(flattmp,2)]);
        QC(i).EPI.SD_demean{j,k}=nanstd(double(tmp2(QC(i).INBRAINMASK.ero0_EPI(:),:)),[],1);
        QC(i).EPI.MEAN_demean{j,k}=nanmean(double(tmp2(QC(i).INBRAINMASK.ero0_EPI(:),:)),1);
        QC(i).EPI.DV{j,k}=rms(diff(flattmp(QC(i).INBRAINMASK.ero0_EPI(:),:),1,2));
        QC(i).EPI.DV{j,k}=[0 QC(i).EPI.DV{j,k}];
    end
    
    % load the motion parameters
    QC(i).motlist{j}= {...
        [QC(i).subdatadir '/pb01.' QC(i).sub '.rest' num2str(j) '.tcat.motion_REF-EPI1.1D'],...
        [QC(i).subdatadir '/pb02.' QC(i).sub '.rest' num2str(j) '.despike.motion_REF-EPI1.1D'],...
        [QC(i).subdatadir '/pb03.' QC(i).sub '.rest' num2str(j) '.despike.tshift.motion_REF-EPI1.1D'],...
        [QC(i).subdatadir '/pb04.' QC(i).sub '.rest' num2str(j) '.tshift.motion_REF-EPI1.1D'],...
        [QC(i).subdatadir '/pb05.' QC(i).sub '.rest' num2str(j) '.tshift.despike.motion_REF-EPI1.1D'],...
        };
    for k=1:numel(QC(i).motlist{j})
        
        % obtain MOTION parameters
        QC(i).MOT.rawAFNI{j,k}=textread(QC(i).motlist{j}{k});
        
        % convert to XYZPRY ordering
        QC(i).MOT.AFNI{j,k}=convert2XYZPRY(QC(i).MOT.rawAFNI{j,k},'afni');
        
        % calculate mean-removed MOTION
        QC(i).MOT.AFNI_toVOL1{j,k}=QC(i).MOT.AFNI{j,k}-repmat((QC(i).MOT.AFNI{j,k}(1,:)),[size(QC(i).MOT.AFNI{j,k},1) 1]);
        
        % calculate FD
        QC(i).FD.AFNI{j,k}=xyzPRY2FD(QC(i).MOT.AFNI{j,k},'degrees');
        
        % calculate ABSTRAN and ABSROT
        QC(i).absTRAN_toVOL1{j,k}=sum(abs(QC(i).MOT.AFNI_toVOL1{j,k}(:,1:3)),2);
        QC(i).absROT_toVOL1{j,k}=sum(abs(QC(i).MOT.AFNI_toVOL1{j,k}(:,4:6)),2);
        
    end
    
     % load the motion parameters
    QC(i).motlistINT{j}= {...
        [QC(i).subdatadir '/pb01.' QC(i).sub '.rest' num2str(j) '.tcat.motion_REF-INT.1D'],...
        [QC(i).subdatadir '/pb02.' QC(i).sub '.rest' num2str(j) '.despike.motion_REF-INT.1D'],...
        [QC(i).subdatadir '/pb03.' QC(i).sub '.rest' num2str(j) '.despike.tshift.motion_REF-INT.1D'],...
        [QC(i).subdatadir '/pb04.' QC(i).sub '.rest' num2str(j) '.tshift.motion_REF-INT.1D'],...
        [QC(i).subdatadir '/pb05.' QC(i).sub '.rest' num2str(j) '.tshift.despike.motion_REF-INT.1D'],...
        };
    for k=1:numel(QC(i).motlistINT{j})
        
        % obtain MOTION parameters
        QC(i).MOT.rawAFNIINT{j,k}=textread(QC(i).motlistINT{j}{k});
        
        % convert to XYZPRY ordering
        QC(i).MOT.AFNIINT{j,k}=convert2XYZPRY(QC(i).MOT.rawAFNIINT{j,k},'afni');
        
        % calculate mean-removed MOTION
        QC(i).MOT.AFNI_toVOL1INT{j,k}=QC(i).MOT.AFNIINT{j,k}-repmat((QC(i).MOT.AFNIINT{j,k}(1,:)),[size(QC(i).MOT.AFNIINT{j,k},1) 1]);
        
        % calculate FD
        QC(i).FD.AFNIINT{j,k}=xyzPRY2FD(QC(i).MOT.AFNIINT{j,k},'degrees');
        
        % calculate ABSTRAN and ABSROT
        QC(i).absTRAN_toVOL1INT{j,k}=sum(abs(QC(i).MOT.AFNI_toVOL1INT{j,k}(:,1:3)),2);
        QC(i).absROT_toVOL1INT{j,k}=sum(abs(QC(i).MOT.AFNI_toVOL1INT{j,k}(:,4:6)),2);
        
    end
    
    % load information about 3dDespike
    QC(i).spikelist{j}= {...
        [QC(i).subdatadir '/pb02.' QC(i).sub '.rest' num2str(j) '.despike_SPIKES.nii.gz'],...
        [QC(i).subdatadir '/pb05.' QC(i).sub '.rest' num2str(j) '.tshift.despike_SPIKES.nii.gz'],...
        };
    QC(i).spikelistbase{j}= {...
        [QC(i).subdatadir '/pb01.' QC(i).sub '.rest' num2str(j) '.tcat4.noskull.nii.gz'],...
        [QC(i).subdatadir '/pb01.' QC(i).sub '.rest' num2str(j) '.tcat4.noskull.nii.gz'],...
        };
    for k=1:numel(QC(i).spikelist{j})
        
        % what's 3dDespike doing? spikes >2.5 are despiked by default
        tmp1=nii_give_RAS(QC(i).spikelist{j}{k}); tmp1d=size(tmp1);
        tmp1=reshape(tmp1,[tmp1d(1)*tmp1d(2)*tmp1d(3) tmp1d(4)]);
        tmp1=tmp1>2.5;
        
        tmp2=nii_give_RAS(QC(i).spikelistbase{j}{k});
        tmp2=~~tmp2;
        QC(i).despike{j,k}=sum(tmp1(tmp2(:),:),1)/nnz(tmp2);
        
    end
    
end

QC2=QC(i);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out264] = give264tcs(flattmp,power264)

[P dims]=nii_give_RAS(power264);

for i=1:size(P,4)
    tmpp=P(:,:,:,i);
    out264(i,:)=nanmean(flattmp(~~tmpp(:),:),1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mot2] = convert2XYZPRY(mot,progtype)

% convert to RAS motion parameters

switch progtype
    
    case 'afni'
        % -YPRxy-z
        mot2=zeros(size(mot));
        mot2(:,1)=mot(:,4);
        mot2(:,2)=mot(:,5);
        mot2(:,3)=-mot(:,6);
        mot2(:,4)=mot(:,2);
        mot2(:,5)=mot(:,3);
        mot2(:,6)=-mot(:,1);
        
    case 'fsl'
        % PRYxyz
        mot2=zeros(size(mot));
        mot2(:,1)=mot(:,4);
        mot2(:,2)=mot(:,5);
        mot2(:,3)=mot(:,6);
        mot2(:,4)=mot(:,1);
        mot2(:,5)=-mot(:,2);
        mot2(:,6)=-mot(:,3);
        
    case 'spm'
        % x-yzP-RY
        mot2=zeros(size(mot));
        mot2(:,1)=mot(:,1);
        mot2(:,2)=-mot(:,2);
        mot2(:,3)=mot(:,3);
        mot2(:,4)=mot(:,4);
        mot2(:,5)=-mot(:,5);
        mot2(:,6)=mot(:,6);
        
    case '4dfp'
        % x-yzP-RY
        mot2=zeros(size(mot));
        mot2(:,1)=-mot(:,1);
        mot2(:,2)=-mot(:,2);
        mot2(:,3)=mot(:,3);
        mot2(:,4)=-mot(:,4);
        mot2(:,5)=mot(:,5);
        mot2(:,6)=mot(:,6);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fd]= xyzPRY2FD(mat,rotkind)

% presume rotations are in degrees

% convert to arc length at 50 mm for rotation

circ=2*pi*50;
switch rotkind
    case 'degrees'
        mat(:,4:6)=circ*mat(:,4:6)/360;
    case 'radians'
        mat(:,4:6)=circ*mat(:,4:6)/(2*pi);
end

fd=sum(abs(diff(mat,1)),2);
fd=[0;fd]; %



