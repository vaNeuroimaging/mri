function JDP_P1_gather(varargin)

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
figdir=[pdir '/FIGS'];
if ~exist(figdir)
    mkdir(figdir);
end
mastermat=[figdir '/P1.mat'];

% where we store summary logs and pictures
[jk logstem jk2]=fileparts(scriptstr);
logdir=['logs_' logstem ];
logpath=[scriptdir '/' logdir ];
clear jk jk2;

% the subjects that could be processed
subjlistfile=[scriptdir '/cleansublist.txt'];
subjlist=textread(subjlistfile,'%s');

ifac=10000;
showfac=100;

[gen]=givegeneralvars();

switch analtype
    
    case 'gather'
        
        isub=[1:numel(subjlist)];
        isub=startsub:endsub;
        
        % gather the QC mats
        for i=isub
            fprintf('Loading P1.mat for subject %d of %d:\t%s\n',i,numel(subjlist),subjlist{i});
            subdatadir=[datadir '/' subjlist{i}];
            load([subdatadir '/P1.mat']);
            QC2(i)=QC;
            clear QC;
        end
        QC=QC2;
        clear QC2;
        
        % read the file that says which runs people will use
        dataselectfile=[logpath '/JDP_P1_preprocess_DATASELECT.txt'];
        [selsub selrun seldecid]=textread(dataselectfile,'%s%d%d');
        ctr=0;
        for i=isub
            cursub=QC(i).sub;
            curruns=QC(i).iruns;
            for k=1:numel(curruns)
                ctr=ctr+1;
                fprintf('Sub %s\t\tRun %d\t\tSelsub %s\tSelrun %d\t Selected %d\n',cursub,curruns(k),selsub{ctr},selrun(ctr),seldecid(ctr));
                if (isequal(cursub,selsub{ctr}) & isequal(curruns(k),selrun(ctr)))
                    QC(i).select.run(k)=seldecid(ctr);
                else
                    fprintf('Mismatch!\n');
                end
                
            end
        end
        
        % determine whether entire subjects not selected
        for i=isub
            if nnz(QC(i).select.run)==0
                QC(i).select.sub=0;
            else
                QC(i).select.sub=1;
            end
            fprintf('Subject %d of %d: %s is used %d\n',i,numel(subjlist),subjlist{i},QC(i).select.sub);
        end
        
        % set up some basic information about run sizes
        for i=isub
            curruns=QC(i).iruns;
            ctr=0;
            for k=1:numel(curruns)
                QC(i).runsize(k)=numel(QC(i).FD.AFNI{1,k});
                QC(i).runborders(k,1)=ctr+1;
                QC(i).runborders(k,2)=QC(i).runborders(k,1)-1+QC(i).runsize(k);
                ctr=QC(i).runborders(k,2);
            end
        end
        
        % set up some basic temporal masks
        
        for i=isub
            curruns=QC(i).iruns;
            for k=1:numel(curruns)
                
                
                
                % denote the first X volumes of each run
                skipvols=str2num(QC(i).SKIPTRS{k}{1}); % skip this many TRs at start of run, will be 8-10 sec
                QC(i).select.runstart{k}=ones(QC(i).runsize(k),1);
                if skipvols>0
                    QC(i).select.runstart{k}(1:skipvols)=0;
                end
                QC(i).select.runstart{k}=~~QC(i).select.runstart{k};
                
                % denote whether selected or not
                if QC(i).select.run(k)==1
                    QC(i).select.runselect{k}=ones(QC(i).runsize(k),1);
                else
                    QC(i).select.runselect{k}=zeros(QC(i).runsize(k),1);
                end
                QC(i).select.runselect{k}=~~QC(i).select.runselect{k};
                
                
                % denote the amount of motion
                QC(i).select.FD050{k}=QC(i).FD.AFNI{1,k}<0.5;
                QC(i).select.FD040{k}=QC(i).FD.AFNI{1,k}<0.4;
                QC(i).select.FD030{k}=QC(i).FD.AFNI{1,k}<0.3;
                QC(i).select.FD020{k}=QC(i).FD.AFNI{1,k}<0.2;
                QC(i).select.FD025{k}=QC(i).FD.AFNI{1,k}<0.25;
            end
        end
        
        % set up some basic QC summaries
        for i=isub
            curruns=QC(i).iruns;
            for k=1:numel(curruns)
                msk=QC(i).select.runstart{k};
                QC(i).FDbar(k)=mean(QC(i).FD.AFNI{1,k}(msk));
                QC(i).FDmedian(k)=median(QC(i).FD.AFNI{1,k}(msk));
                QC(i).DVbar(k)=mean(QC(i).EPI.DV{k,1}(msk));
                QC(i).DVmedian(k)=median(QC(i).EPI.DV{k,1}(msk));
                
            end
        end
        
        save(mastermat,'QC','-v7.3');
     
        
    case 'reload_selection'
        
        QC=varargin{1,6};
        dataselectfile=varargin{1,7};
        
        [selsub selrun seldecid]=textread(dataselectfile,'%s%d%d');
        ctr=0;
        for i=1:numel(QC)
            cursub=QC(i).sub;
            curruns=QC(i).iruns;
            for k=1:numel(curruns)
                ctr=ctr+1;
                fprintf('Sub %s\t\tRun %d\t\tSelsub %s\tSelrun %d\t Selected %d\n',cursub,curruns(k),selsub{ctr},selrun(ctr),seldecid(ctr));
                if (isequal(cursub,selsub{ctr}) & isequal(curruns(k),selrun(ctr)))
                    QC(i).select.run(k)=seldecid(ctr);
                else
                    fprintf('Mismatch!\n');
                end
                
            end
        end
       
       save(mastermat,'QC','-v7.3');
        
       
    case 'showtotaltimes'
        
        QC=varargin{1,6};
        [QC] = addvariables(QC);
        
        odir = [figdir '/NP1_F2'];
        mkdir(odir);
        
        % constrain each plot to the same length of time
        % i.e., shows constant time scale
        tlimited=0;
        tshown=16; % minutes at a time on a plot (if tlimited)
        
        isub=[1:numel(subjlist)];
        isub=startsub:endsub;
        
        for i=isub
            
%             fprintf('Sub %d\t %s\n',i,QC(i).sub);
            catvec=find(QC(i).select.run);
            if ~isempty(catvec)
                
                tmp=catQC(QC,i,catvec);
                
                fprintf('Sub %d\t %s\t\t%2.2f\n',i,QC(i).sub,tmp.totaltime);
                
            end
        end
             
        
    case 'NP1_F2_QCtraces_allruns'
        
        QC=varargin{1,6};
        
        [QC] = addvariables(QC);
        
        odir = [figdir '/NP1_F2'];
        if ~exist(odir)
            mkdir(odir);
        end
        
        % constrain each plot to the same length of time
        % i.e., shows constant time scale
        tlimited=1;
        tshown=17; % minutes at a time on a plot (if tlimited)
        
        % do all runs or just selected?
        doall=1;
        
        isub=startsub:endsub;
        
        for i=isub
            if doall
                catvec=QC(i).iruns;
            else
                catvec=find(QC(i).select.run);
            end
            if ~isempty(catvec)
                
                for k=catvec
                    tmask{k}=~QC(i).select.runstart{k};
%                     tmask{k}=zeros(1,numel(QC(i).FD.AFNI{k,1}));
%                     skptrs=(QC(i).SKIPTRS{k}{1});
%                     if skptrs
%                         tmask{k}(1:skptrs)=1;
%                     end
                end
                
                tmp=catQC(QC,i,catvec,1,tmask);
                tottmask=cat(1,tmask{catvec});
                
                tottmask=~tottmask;
                % calculate some useful numbers
                tmpc=corrcoef([tmp.EPI.GM_mean(tottmask)' tmp.EPI.WM_mean(tottmask)']);
                corrs.GMWM=tmpc(1,2);
                tmpc=corrcoef([tmp.EPI.GM_mean(tottmask)' tmp.EPI.CSF_mean(tottmask)']);
                corrs.GMCSF=tmpc(1,2);
                corrs.rmsmot=rms(rms(tmp.MOT.AFNI(tottmask,:)));
                corrs.meanFD=(mean(tmp.FD.AFNI(tottmask)));
                corrs.medFD=(median(tmp.FD.AFNI(tottmask)));
                corrs.FDover2=100*(nnz(tmp.FD.AFNI(tottmask)>.3)/nnz(tottmask));
                corrs.DSPKover10=100*(nnz(tmp.despike(tottmask)>.10)/nnz(tottmask));
                
                % put nan in unused TRs to show censoring/skipped TRs
                tmp.MOT.AFNI(~tottmask,:)=nan;
                tmp.FD.AFNI(~tottmask)=nan;
                tmp.absTRAN_toVOL1(~tottmask)=nan;
                tmp.absROT_toVOL1(~tottmask)=nan;
                tmp.despike(~tottmask)=nan;
                tmp.EPI.DV(~tottmask)=nan;
                tmp.EPI.CSF_mean(~tottmask)=nan;
                tmp.EPI.WM_mean(~tottmask)=nan;
                tmp.EPI.GM_mean(~tottmask)=nan;
              
                
                fprintf('Sub %d\t %s\t\t%2.2f\n',i,QC(i).sub,tmp.totaltime);

                % assumes all runs have same TR
                trper=ceil(tshown*60/tmp.runborders(1,4));
                if ~tlimited
                    trper=tmp.totaltrs; % uncomment this stop tshown
                end
                tshown2=trper*tmp.runborders(1,3)/60;
                vborders=[tmp.runborders(2:end,1); tmp.runborders(end,2)];
                
                showset=1:trper:tmp.totaltrs;
                if numel(showset)==1
                    showxlim=[1 trper];
                    showxplt=[1 tmp.totaltrs];
                else
                    showset=[showset showset(end)+trper];
                    showxlim=[[showset(1:end-1)]' [showset(2:end)-1]'];
                    showxplt=showxlim;
                    showxplt(end,end)=tmp.totaltrs;
                end
                
                xoff=.01;
                xpos=[0.9 0.8 0.7]+xoff;
                ypos=0.007;
                
                for shows=1:size(showxlim,1)
                    sxl=[showxplt(shows,1):showxplt(shows,2)];
                    sxlimz=showxlim(shows,:);
%                     tvborders=vborders(ismember(vborders,ax1));
                    
                    
                    close all;
                    h=figure;
                    set(h,'position',[10 10 1920 1080]);
                    set(h,'visible','off');
                    tfs=16;
                    lfs=10;
                    ylpos=-0.025;
                    lw=2;
                    
                    subplot(5,1,1);
                    ylimz=[-1 1];
                    moffset=[.25 .15 .05 -.05 -.15 -.25];
                    hold on;
                    for j=1:6
%                         plot(sxl,tmp.MOT.AFNI(sxl,j)+moffset(j),'color',gen.mclrs(j,:),'linewidth',lw);
                        plot(sxl,tmp.MOT.AFNI(sxl,j),'color',gen.mclrs(j,:),'linewidth',lw);
                    end
                    hold off
                    ylim(ylimz);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',ylimz,'box','on','fontsize',tfs);
                    hh=legend({'x','y','z','P','R','Y'},'location','NorthEast');
                    set(hh,'fontsize',lfs);
                    title({['Subject ' tmp.sub ],['Run: ' num2str(catvec)]},'Interpreter','none');
                    yl=ylabel('mm');
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    vline(vborders,'k');
                    str=['RMS motion: ' num2str(corrs.rmsmot,'%3.3f') ];
                    text(ypos,xpos(1),str,'units','normalized','horizontalalignment','left');
                    
                    subplot(5,1,2);
                    ylimz=[0 2];
                    hold on;
                    plot(sxl,tmp.FD.AFNI(sxl),'color',gen.fdclr,'linewidth',lw);
                    plot(sxl,tmp.absTRAN_toVOL1(sxl),':','color',gen.abstranclr,'linewidth',lw);
                    plot(sxl,tmp.absROT_toVOL1(sxl),'--','color',gen.absrotclr,'linewidth',lw)
                    hold off;
                    ylim(ylimz);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',ylimz,'box','on','fontsize',tfs);
                    hh=legend({'FD','absT','absR'},'location','NorthEast');
                    set(hh,'fontsize',lfs);
                    yl=ylabel('mm');
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    vline(vborders,'k');
                    str=['mean FD: ' num2str(corrs.meanFD,'%3.3f') ];
                    text(ypos,xpos(1),str,'units','normalized','horizontalalignment','left');
                    str=['med   FD: ' num2str(corrs.medFD,'%3.3f') ];
                    text(ypos,xpos(2),str,'units','normalized','horizontalalignment','left');
                    str=['%>0.30 mm: ' num2str(corrs.FDover2,'%3.1f') ];
                    text(ypos,xpos(3),str,'units','normalized','horizontalalignment','left');
                    
                    subplot(5,1,3);
                    ylimz=[0 50];
                    hold on;
                    plot(sxl,tmp.EPI.DV(sxl),'color',gen.dvclr,'linewidth',lw);
                    hold off;
                    ylim(ylimz);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',ylimz,'box','on','fontsize',tfs);
                    hh=legend({'DV'},'location','NorthEast');
                    set(hh,'fontsize',lfs);
                    yl=ylabel('signal (a.u.)');
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    vline(vborders,'k');
                    
                    subplot(5,1,4);
                    ylimz=[0 100];
                    hold on;
                    plot(sxl,tmp.despike(sxl)*100,'color',gen.dspkclr,'linewidth',lw);
                    hold off;
                    ylim(ylimz);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',ylimz,'box','on','fontsize',tfs);
                    hh=legend({'3dDespike'},'location','NorthEast');
                    set(hh,'fontsize',lfs);
                    yl=ylabel('% image despiked');
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    vline(vborders,'k');
                    str=['%>10%: ' num2str(corrs.DSPKover10,'%3.1f') ];
                    text(ypos,xpos(1),str,'units','normalized','horizontalalignment','left');
                    
                    subplot(5,1,5);
                    limz=[-20 20];
                   
                    imagesc(tmp.EPI.Power264(:,:),limz);
                    colormap(gray);
                    xlim(sxlimz);
                    set(gca,'xtick',sxlimz,'ytick',[],'box','on','fontsize',tfs);
                    xlabel({['Volume #'],['(TR=' num2str(tmp.runborders(1,4),'%2.2f') ' sec ; total=' num2str(tmp.totaltime,'%2.2f') ' min; ' num2str(tshown,'%2.2f') ' min shown)']});
                    yl=ylabel('264 ROIs');
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    
                    
                    set(h,'paperpositionmode','auto');
                    ofile=[odir '/' num2str(i*ifac+shows*showfac+0) '.png'];
                    print(gcf,'-dpng',ofile);
                    close(h);
                end
                
                for shows=1:size(showxlim,1)
                    sxl=[showxplt(shows,1):showxplt(shows,2)];
                    sxlimz=showxlim(shows,:);
                    
                    
                    close all;
                    h=figure;
                    set(h,'position',[10 10 1920 1080]);
                    set(h,'visible','off');
                    tfs=16;
                    lfs=10;
                    ylpos=-0.025;
                    lw=2;
                    
                    
                    subplot(5,1,1);
                    ylimz=[-1 1];
                    hold on;
                    for j=1:6
                        plot(sxl,tmp.MOT.AFNI(sxl,j),'color',gen.mclrs(j,:),'linewidth',lw);
                    end
                    hold off
                    ylim(ylimz);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',ylimz,'box','on','fontsize',tfs);
                    hh=legend({'x','y','z','P','R','Y'},'location','NorthEast');
                    set(hh,'fontsize',lfs);
                    title({['Subject ' tmp.sub ],['Run: ' num2str(catvec)]},'Interpreter','none');
                    yl=ylabel('mm');
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    vline(vborders,'k');
                    str=['RMS motion: ' num2str(corrs.rmsmot,'%3.3f') ];
                    text(ypos,xpos(1),str,'units','normalized','horizontalalignment','left');
                    
                    subplot(5,1,2);
                    gmoff=10;
                    wmoff=0;
                    csfoff=-10;
                    ylimz=[-30 30];
                    hold on;
                   
                   
                    plot(sxl,tmp.EPI.GM_mean(sxl)+gmoff,'color',gen.meanclr,'linewidth',lw);
                     plot(sxl,tmp.EPI.WM_mean(sxl)+wmoff,'color',gen.wmclr,'linewidth',lw);
                      plot(sxl,tmp.EPI.CSF_mean(sxl)+csfoff,'color',gen.csfclr,'linewidth',.5);
%                     plot(sxl,tmp.EPI.GM_std(sxl),'color',gen.stdclr,'linewidth',lw);
                    hold off
                    ylim(ylimz);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',ylimz,'box','on','fontsize',tfs);
                    hh=legend({'GM','WM','CSF'},'location','NorthEast');
                    set(hh,'fontsize',lfs);
                    yl=ylabel('signal (a.u.)');
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    vline(vborders,'k');
                     str=['GM:WM  r: ' num2str(corrs.GMWM,'%3.2f') ];
                    text(ypos,xpos(1),str,'units','normalized','horizontalalignment','left');
                    str=['GM:CSF r: ' num2str(corrs.GMCSF,'%3.2f') ];
                    text(ypos,xpos(2),str,'units','normalized','horizontalalignment','left');
                    
                    
                    subplot(5,1,3);
                    limz=[-20 20];
                    imagesc(tmp.EPI.CSF,limz);
                    colormap(gray);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',[],'box','on','fontsize',tfs);
                    yl=ylabel(['CSF voxels (' num2str(size(tmp.EPI.CSF,1)) ')']);
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    
                    subplot(5,1,4);
                    limz=[-20 20];
                    imagesc(tmp.EPI.WM,limz);
                    colormap(gray);
                    xlim(sxlimz);
                    set(gca,'xtick',[],'ytick',[],'box','on','fontsize',tfs);
                    yl=ylabel(['WM voxels (' num2str(size(tmp.EPI.WM,1)) ')']);
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    
                    subplot(5,1,5);
                    limz=[-20 20];
                    imagesc(tmp.EPI.GM,limz);
                    colormap(gray);
                    xlim(sxlimz);
                    set(gca,'xtick',sxlimz,'ytick',[],'box','on','fontsize',tfs);
                    xlabel({['Volume #'],['(TR=' num2str(tmp.runborders(1,4),'%2.2f') ' sec ; total=' num2str(tmp.totaltime,'%2.2f') ' min; ' num2str(tshown,'%2.2f') ' min shown)']});
                    yl=ylabel(['GM voxels (' num2str(size(tmp.EPI.GM,1)) ')']);
                    set(yl,'units','normalized','position',[ylpos .5 0]);
                    
                    set(h,'paperpositionmode','auto');
                    ofile=[odir '/' num2str(i*ifac+shows*showfac+1) '.png'];
                    print(gcf,'-dpng',ofile);
                    
                    close(h);
                    
                end
                
                
                
            end
            
        end
        
        
        
end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gen] = givegeneralvars()

gen.mclrs=[
    255 0 0; %X - red
    0 255 0; %Y - green
    0 0 255; %Z - blue
    0 255 255; %Pitch - light blue
    255 0 255; %Roll - pink
    220 220 0 %Yaw - yellow
    ]/255;

gen.abstranclr=[.5 0 0];
gen.absrotclr=[.5 0 0];
gen.dvclr=[0 0 .5];
gen.dspkclr=[.5 .5 .5];
gen.fdclr=[1 0 0];
gen.stdclr=[0 .5 0];
gen.meanclr=[0 0 0];
gen.wmclr=[1 .5 0];
gen.csfclr=[0.5 0.25 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [brd] = giveborders(szvec)

runends=cumsum(szvec);
runstarts(1)=1; runstarts(2:numel(runends))=runends(1:end-1)+1;

brd=[runstarts' runends'];

brd(:,3)=brd(:,2)-brd(:,1)+1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmp] = catQC(QC,i,catvec,ord,tmask)

% concatenates across runs
% also applies mode 1000 to EPI prior to concatenatinos

% remove mean and linear terms

for k=catvec
    skiptr=str2num(QC(i).SKIPTRS{k}{1})+1;
    mtc{k}=QC(i).MOT.AFNI_toVOL1{k,1}-repmat(QC(i).MOT.AFNI_toVOL1{k,1}(skiptr+1,:),[size(QC(i).MOT.AFNI_toVOL1{k,1},1) 1]);
    mtc2{k}=sum(abs(mtc{k}(:,1:3)),2);
    mtc3{k}=sum(abs(mtc{k}(:,4:6)),2);
end
tmp.MOT.AFNI=cat(1,mtc{catvec});
tmp.absTRAN_toVOL1=cat(1,mtc2{catvec});
tmp.absROT_toVOL1=cat(1,mtc3{catvec});
% tmp.MOT.AFNI=cat(1,QC(i).MOT.AFNI_toVOL1{catvec,1});
% tmp.absTRAN_toVOL1=cat(1,QC(i).absTRAN_toVOL1{catvec,1});
% tmp.absROT_toVOL1=cat(1,QC(i).absROT_toVOL1{catvec,1});

tmp.FD.AFNI=cat(1,QC(i).FD.AFNI{catvec,1});
tmp.despike=[cat(2,QC(i).despike{catvec,1})];

for k=catvec
    QC(i).EPI.DV{k,1}=QC(i).EPI.DV{k,1}*QC(i).mode1000{k};
end
tmp.EPI.DV=[cat(2,QC(i).EPI.DV{catvec,1})];


for k=catvec
    QC(i).EPI.Power264{k,1}=QC(i).EPI.Power264{k,1}*QC(i).mode1000{k};
    mtc{k}=removepoly(QC(i).EPI.Power264{k,1},ord,tmask{k});
end
tmp.EPI.Power264=cat(2,mtc{catvec});

for k=catvec
    QC(i).EPI.CSF{k,1}=QC(i).EPI.CSF{k,1}*QC(i).mode1000{k};
    mtc{k}=removepoly(QC(i).EPI.CSF{k,1},ord,tmask{k});
    mtc_mean{k}=mean(mtc{k},1);
end
tmp.EPI.CSF=cat(2,mtc{catvec});
tmp.EPI.CSF_mean=cat(2,mtc_mean{catvec});

for k=catvec
    QC(i).EPI.WM{k,1}=QC(i).EPI.WM{k,1}*QC(i).mode1000{k};
    mtc{k}=removepoly(QC(i).EPI.WM{k,1},ord,tmask{k});
    mtc_mean{k}=mean(mtc{k},1);
end
tmp.EPI.WM=cat(2,mtc{catvec});
tmp.EPI.WM_mean=cat(2,mtc_mean{catvec});

for k=catvec
    QC(i).EPI.GM{k,1}=QC(i).EPI.GM{k,1}*QC(i).mode1000{k};
    mtc{k}=removepoly(QC(i).EPI.GM{k,1},ord,tmask{k});
    mtc_mean{k}=mean(mtc{k},1);
    mtc_std{k}=std(mtc{k},[],1);
end
tmp.EPI.GM=cat(2,mtc{catvec});
tmp.EPI.GM_mean=cat(2,mtc_mean{catvec});
tmp.EPI.GM_std=cat(2,mtc_std{catvec});

for k=catvec
    QC(i).EPI.MEAN_demean{k,1}=QC(i).EPI.MEAN_demean{k,1}*QC(i).mode1000{k};
end
tmp.EPI.MEAN_demean=[cat(2,QC(i).EPI.MEAN_demean{catvec,1})];

for k=catvec
    QC(i).EPI.SD_demean{k,1}=QC(i).EPI.SD_demean{k,1}*QC(i).mode1000{k};
end
tmp.EPI.SD_demean=[cat(2,QC(i).EPI.SD_demean{catvec,1})];

tmp.sub=QC(i).sub;

[tmp.runborders]=giveborders(QC(i).runsize(catvec));
for ctr=1:numel(catvec)
    tmp.runborders(ctr,4)=QC(i).TR{catvec(ctr)};
    tmp.runborders(ctr,5)=QC(i).TR{catvec(ctr)}*tmp.runborders(ctr,3);
    tmp.runborders(ctr,6)=tmp.runborders(ctr,5)/60;
end
tmp.totaltime=sum(tmp.runborders(:,6));
tmp.totaltrs=sum(tmp.runborders(:,3));







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [QC] = addvariables(QC)


% ______________
% **************
% ^^^^^^^^^^^^^^
% mode 1000 values; many possible ways to calculate this

isub=[1:numel(QC)];
for i=isub
    for k=QC(i).iruns
        tmp=QC(i).EPI.ANATAVE{k,1};
        tmp2=QC(i).GMMASK.RIBBONMASK_ero0_EPI;
        gmvals=tmp(tmp2);
        gmvals(gmvals<10)=[];
        [h hh]=hist(gmvals,100);
        QC(i).mode1000{k}=1000./(hh(find(h==max(h))));
        if numel(QC(i).mode1000{k})>1
            QC(i).mode1000{k}=mean(QC(i).mode1000{k});
        end

    end
end

% ______________
% **************
% ^^^^^^^^^^^^^^
% zero POSITION variables to TR5; was a bad choice to do it to 1
% fine for WU data but not any other datasets. 
% 
% isub=[1:numel(QC)];
% for i=isub
%     for k=QC(i).iruns
%         tmp=QC(i).EPI.ANATAVE{k,1};
%         tmp2=QC(i).GMMASK.RIBBONMASK_ero0_EPI;
%         gmvals=tmp(tmp2);
%         gmvals(gmvals<10)=[];
%         [h hh]=hist(gmvals,100);
%         QC(i).mode1000{k}=1000./(hh(find(h==max(h))));
%         if numel(QC(i).mode1000{k})>1
%             QC(i).mode1000{k}=mean(QC(i).mode1000{k});
%         end
% 
%     end
% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [resid pred b] = removepoly(tc,ord,varargin)

% remove legendre polynomials up to ORD
% 0 - mean
% 1 - trend
% 2 - bow down 
% 3 - sinusoid
% 4 - w shape
% 5 - w sinusoid
% etc etc
%
% be careful, these are high-pass filters also

% presumes tc is vox x time

d=size(tc);

for i=0:ord
    b=legendre(i,linspace(-1,1,d(2)));
    if i==0
        r=b(1,:);
    else
        r=[r; b(1,:)];
    end
end


if isempty(varargin)
    % use all timepoints
    b=r'\tc';
    pred=r'*b;
    resid=tc-pred';
    
else
    % use only specified timepoints
    tmask=varargin{1,1};
    b=r(:,~tmask)'\tc(:,~tmask)';
    pred=r'*b;
    resid=tc-pred';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [clr2]=lighten(clr,gain)

d=size(clr);
clr2=(ones(d)-clr)*gain+clr;

    
