
threshold = 5.0;

data = ['/fmri/data3/Evan/Gene-Rest-Nback//Analysis/ICA/Old/ICA_lowthresh20dim_HBMpaper.gica/groupmelodic.ica/stats/thresh_zstat3.img'];


load /fmri/data3/Evan/Gene-Rest-Nback/Scripts/FSLpeaks_template.mat


%outputname = [outputpath subj '_' filetypes{filetype} '_roi.mat'];

o = [];
d = [];

imgname = data;

[p f e] = fileparts(imgname);
binf = 1;

func = ['img >= ' num2str(threshold)];

d = f; l = f;
if ~isempty(func)
    d = [d ' func: ' func];
    l = [l '_f_' func];
end
if binf
    d = [d ' - binarized'];
    l = [l '_bin'];
end
o = maroi_image(struct('vol', spm_vol(imgname), 'binarize',binf,...
    'func', func));

% convert to matrix format to avoid delicacies of image format
o = maroi_matrix(o);


o = descrip(o,d);
o = label(o,l);

%varargout = {saveroi(o, outputname)};

[Y multv vXYZ mat] = getdata(o, data, 'l');

xSPM.Z = Y;
xSPM.XYZ = vXYZ;
xSPM.M = mat;
xSPM.XYZmm = mat(1:3,1:3)*vXYZ + repmat(mat(1:3,4),1,size(vXYZ,2));
xSPM.STAT = 'Z';
xSPM.STATstr = 'Z';
xSPM.title = 'Image';
xSPM.u = threshold;
xSPM.thresDesc = 'User-set';
xSPM.k = 0;


%[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');

%hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
%hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,xSPM.M,[44;53;44],1);


%spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');


%-Initialise
%----------------------------------------------------------------------
%SPMid      = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Results');
FS         = spm('FontSizes');

% clear satfig if it exists
%----------------------------------------------------------------------
hSat       = findobj('tag','Satellite');
spm_figure('clear',hSat);

%-GexSPM.STATt thresholded xSPM data and parameters of design
%=======================================================================
%     if (nargin > 1)
% 	[SPM,xSPM] = spm_getSPM(varargin{2});
%     else
% 	[SPM,xSPM] = spm_getSPM;
%     end;

if isempty(xSPM)
    varargout = {[],[],[]};
    return;
end

M         = xSPM.M;
DIM       = xSPM.DIM;
try
    units = xSPM.units;
catch
    units = {'mm' 'mm' 'mm'};
end

% ensure pwd = swd so that relative filenames are valid
%----------------------------------------------------------------------
%cd(SPM.swd)

%-Setup Results User Interface; Display MIP, design matrix & parameters
%======================================================================
spm('FigName',['SPM{',xSPM.STAT,'}: Results'],Finter,CmdLine);


%-Setup results GUI
%----------------------------------------------------------------------
spm_figure('Clear',Finter)
hReg      = spm_results_ui('SetupGUI',M,DIM,xSPM,Finter);

%-Setup design interrogation menu
%----------------------------------------------------------------------
%hDesRepUI = spm_DesRep('DesRepUI',SPM);
figure(Finter)

%-Setup Maximium intensity projection (MIP) & register
%----------------------------------------------------------------------
hMIPax = axes('Parent',Fgraph,'Position',[0.05 0.60 0.55 0.36],'Visible','off');
hMIPax = spm_mip_ui(xSPM.Z,xSPM.XYZmm,M,DIM,hMIPax,units);

spm_XYZreg('XReg',hReg,hMIPax,'spm_mip_ui');
if xSPM.STAT == 'P'
    str = xSPM.STATstr;
else
    str = ['SPM\{',xSPM.STATstr,'\}'];
end
text(240,260,str,...
    'Interpreter','TeX',...
    'FontSize',FS(14),'Fontweight','Bold',...
    'Parent',hMIPax)


%-Print comparison title
%----------------------------------------------------------------------
hTitAx = axes('Parent',Fgraph,...
    'Position',[0.02 0.95 0.96 0.02],...
    'Visible','off');

text(0.5,0,xSPM.title,'Parent',hTitAx,...
    'HorizontalAlignment','center',...
    'VerticalAlignment','baseline',...
    'FontWeight','Bold','FontSize',FS(14))


%-Print SPMresults: Results directory & thresholding info
%----------------------------------------------------------------------
hResAx = axes('Parent',Fgraph,...
    'Position',[0.05 0.55 0.45 0.05],...
    'DefaultTextVerticalAlignment','baseline',...
    'DefaultTextFontSize',FS(9),...
    'DefaultTextColor',[1,1,1]*.7,...
    'Units','points',...
    'Visible','off');
AxPos = get(hResAx,'Position'); set(hResAx,'YLim',[0,AxPos(4)])
h     = text(0,24,'SPMresults:','Parent',hResAx,...
    'FontWeight','Bold','FontSize',FS(14));
text(get(h,'Extent')*[0;0;1;0],24,spm_str_manip(data,'a30'),'Parent',hResAx)
try
    thresDesc = xSPM.thresDesc;
    text(0,12,sprintf('Height threshold %c = %0.6f  {%s}',xSPM.STAT,xSPM.u,thresDesc),'Parent',hResAx)
catch
    text(0,12,sprintf('Height threshold %c = %0.6f',xSPM.STAT,xSPM.u),'Parent',hResAx)
end
text(0,00,sprintf('Extent threshold k = %0.0f voxels',xSPM.k), 'Parent',hResAx)


% %-Plot design matrix
% %----------------------------------------------------------------------
% hDesMtx   = axes('Parent',Fgraph,'Position',[0.65 0.55 0.25 0.25]);
% hDesMtxIm = image((SPM.xX.nKX + 1)*32);
% xlabel('Design matrix')
% set(hDesMtxIm,'ButtonDownFcn','spm_DesRep(''SurfDesMtx_CB'')',...
%     'UserData',struct(...
%     'X',		SPM.xX.xKXs.X,...
%     'fnames',	{reshape({SPM.xY.VY.fname},size(SPM.xY.VY))},...
%     'Xnames',	{SPM.xX.name}))
% 
% %-Plot contrasts
% %----------------------------------------------------------------------
% nPar   = size(SPM.xX.X,2);
% xx     = [repmat([0:nPar-1],2,1);repmat([1:nPar],2,1)];
% nCon   = length(xSPM.Ic);
% xCon   = SPM.xCon;
% if nCon
%     dy     = 0.15/max(nCon,2);
%     hConAx = axes('Position',[0.65 (0.80 + dy*.1) 0.25 dy*(nCon-.1)],...
%         'Tag','ConGrphAx','Visible','off');
%     title('contrast(s)')
%     htxt   = get(hConAx,'title');
%     set(htxt,'Visible','on','HandleVisibility','on')
% end
% 
% for ii = nCon:-1:1
%     axes('Position',[0.65 (0.80 + dy*(nCon - ii +.1)) 0.25 dy*.9])
%     if xCon(xSPM.Ic(ii)).STAT == 'T' & size(xCon(xSPM.Ic(ii)).c,2) == 1
% 
%         %-Single vector contrast for SPM{t} - bar
%         %--------------------------------------------------------------
%         yy = [zeros(1,nPar);repmat(xCon(xSPM.Ic(ii)).c',2,1);zeros(1,nPar)];
%         h  = patch(xx,yy,[1,1,1]*.5);
%         set(gca,'Tag','ConGrphAx',...
%             'Box','off','TickDir','out',...
%             'XTick',spm_DesRep('ScanTick',nPar,10) - 0.5,'XTickLabel','',...
%             'XLim',	[0,nPar],...
%             'YTick',[-1,0,+1],'YTickLabel','',...
%             'YLim',[min(xCon(xSPM.Ic(ii)).c),max(xCon(xSPM.Ic(ii)).c)] +...
%             [-1 +1] * max(abs(xCon(xSPM.Ic(ii)).c))/10	)
% 
%     else
% 
%         %-F-contrast - image
%         %--------------------------------------------------------------
%         h = image((xCon(xSPM.Ic(ii)).c'/max(abs(xCon(xSPM.Ic(ii)).c(:)))+1)*32);
%         set(gca,'Tag','ConGrphAx',...
%             'Box','on','TickDir','out',...
%             'XTick',spm_DesRep('ScanTick',nPar,10),'XTickLabel','',...
%             'XLim',	[0,nPar]+0.5,...
%             'YTick',[0:size(SPM.xCon(xSPM.Ic(ii)).c,2)]+0.5,....
%             'YTickLabel','',...
%             'YLim',	[0,size(xCon(xSPM.Ic(ii)).c,2)]+0.5	)
% 
%     end
%     ylabel(num2str(xSPM.Ic(ii)))
%     set(h,'ButtonDownFcn','spm_DesRep(''SurfCon_CB'')',...
%         'UserData',	struct(	'i',		xSPM.Ic(ii),...
%         'h',		htxt,...
%         'xCon',		xCon(xSPM.Ic(ii))))
% endxSPM.k


%-Store handles of results section Graphics window objects
%----------------------------------------------------------------------
H  = get(Fgraph,'Children');
H  = findobj(H,'flat','HandleVisibility','on');
H  = findobj(H);
Hv = get(H,'Visible');
set(hResAx,'Tag','PermRes','UserData',struct('H',H,'Hv',{Hv}))


%-Finished results setup
%----------------------------------------------------------------------
%varargout = {hReg,xSPM,SPM};
spm('Pointer','Arrow')


spm_list('List',xSPM,hReg);