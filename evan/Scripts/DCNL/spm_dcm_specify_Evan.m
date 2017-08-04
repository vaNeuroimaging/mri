function DCM = spm_dcm_specify
% Specify inputs of a DCM
% FORMAT [DCM] = spm_dcm_specify
%
% DCM  - the DCM structure (see spm_dcm_ui)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_dcm_specify.m 2671 2009-01-29 19:37:36Z klaas $

Finter = spm_figure('GetWin','Interactive');
WS     = spm('WinScale');

%==========================================================================
% Get design and directory
%==========================================================================
[spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
if ~sts, DCM = []; return; end
swd = spm_str_manip(spmmatfile,'H');
try
    load(fullfile(swd,'SPM.mat'))
catch
    error(['Cannot read ' fullfile(swd,'SPM.mat')]);
end
    
%==========================================================================
% Name
%==========================================================================
name  = spm_input('name for DCM_???.mat','+1','s');

%==========================================================================
% Outputs
%==========================================================================

%-Get cell array of region structures
%--------------------------------------------------------------------------
P     = cellstr(spm_select([2 8],'^VOI.*\.mat$',{'select VOIs'},'',swd));
m     = numel(P);
for i = 1:m
    p     = load(P{i},'xY','-mat');
    xY(i) = p.xY;
end

%==========================================================================
% Inputs
%==========================================================================

%-Get 'causes' or inputs U
%--------------------------------------------------------------------------
spm_input('Input specification:...  ',1,'d');
Sess   = SPM.Sess(xY(1).Sess);
U.dt   = Sess.U(1).dt;
u      = length(Sess.U);
U.name = {};
U.u    = [];
for  i = 1:u
    for  j = 1:length(Sess.U(i).name)
        str   = ['include ' Sess.U(i).name{j} '?'];
        if spm_input(str,2,'y/n',[1 0])
            U.u             = [U.u Sess.U(i).u(33:end,j)];
            U.name{end + 1} = Sess.U(i).name{j};
        end
    end
end

%==========================================================================
% Graph connections
%==========================================================================
n     = size(U.u,2);
a     = zeros(m,m);
c     = zeros(m,n);
b     = zeros(m,m,n);
q     = uicontrol(Finter,'String','done','Position',[300 100 060 020].*WS);
dx    = 20;

%-Intrinsic connections (A matrix)
%--------------------------------------------------------------------------
spm_input('Specify intrinsic connections from',1,'d')
spm_input('to',3,'d')

for i = 1:m
    str    = sprintf('%s %i',xY(i).name,i);
    h1(i)  = uicontrol(Finter,'String',str,...
        'Style','text',...
        'HorizontalAlignment','right',...
        'Position',[080 356-dx*i 080 020].*WS);
    h2(i)  = uicontrol(Finter,'String',sprintf('%i',i),...
        'Style','text',...
        'Position',[180+dx*i 356 020 020].*WS);
end
for i = 1:m
    for j = 1:m
        h3(i,j) = uicontrol(Finter,...
            'Position',[180+dx*j 360-dx*i 020 020].*WS,...
            'Style','radiobutton');
        if i == j
            set(h3(i,j),'Value',1,...
                'enable','off',...
                'BackgroundColor',[0.5 0.5 0.5]);
        end

    end
end
drawnow

%-Wait for 'done'
%--------------------------------------------------------------------------
while true
    pause(0.01);
    if strcmp(get(gco,'Type'),'uicontrol')
        if strcmp(get(gco,'String'),'done')
            for i = 1:m
                for j = 1:m
                    a(i,j) = get(h3(i,j),'Value');
                end
            end
            delete([h1(:); h2(:); h3(:)]);
            break
        end
    end
end

%-Effects of causes (B and C matrices)
%--------------------------------------------------------------------------
for k = 1:n

    %-Buttons and labels
    %----------------------------------------------------------------------
    str   = sprintf(...
        'Effects of %-12s on regions... and connections',...
        U.name{k});
    spm_input(str,1,'d');

    dx    = 20;
    for i = 1:m
        h1(i)  = uicontrol(Finter,'String',xY(i).name,...
            'Style','text',...
            'Position',[080 356-dx*i 080 020].*WS);
        h2(i)  = uicontrol(Finter,...
            'Position',[160 360-dx*i 020 020].*WS,...
            'Style','radiobutton');
    end
    for i = 1:m
        for j = 1:m
            if a(i,j)==1
                % If there is an intrinsic connection
                % allow it to be modulated
                h3(i,j) = uicontrol(Finter,...
                    'Position',[220+dx*j 360-dx*i 020 020].*WS,...
                    'Style','radiobutton');
            end
        end
    end
    drawnow

    %-Wait for 'done'
    %----------------------------------------------------------------------
    set(gcf,'CurrentObject',h2(1))
    while(1)
        pause(0.01)
        if strcmp(get(gco,'Type'),'uicontrol')
            if strcmp(get(gco,'String'),'done')

                %-Get c
                %----------------------------------------------------------
                for i = 1:m
                    c(i,k)   = get(h2(i)  ,'Value');
                end

                %-Get b allowing any 2nd order effects
                %----------------------------------------------------------
                for i = 1:m
                    for j = 1:m
                        if a(i,j)==1
                            b(i,j,k) = get(h3(i,j),'Value');
                        end
                    end
                end
                delete([h1(:); h2(:); h3(find(a==1))])
                break

            end
        end
    end
end
delete(q)


%-Effects of nonlinear modulations (D matrices)
%--------------------------------------------------------------------------
str   = ['Nonlinear DCM ?'];
if spm_input(str,2,'y/n',[1 0])
    nlDCM = 1;
    q     = uicontrol(Finter,'String','done','Position',[300 100 060 020].*WS);
    for k = 1:m

        %-Buttons and labels
        %----------------------------------------------------------------------
        str   = sprintf(...
            'Effects of %-12s activity on connections',...
            xY(k).name);
        spm_input(str,1,'d');

        dx    = 20;
        for i = 1:m
            for j = 1:m
                if a(i,j)==1
                    % If there is an intrinsic connection
                    % allow it to be modulated
                    h4(i,j) = uicontrol(Finter,...
                        'Position',[220+dx*j 360-dx*i 020 020].*WS,...
                        'Style','radiobutton');
                end
            end
        end
        drawnow

        %-Wait for 'done'
        %----------------------------------------------------------------------
        set(gcf,'CurrentObject',h4(1))
        while(1)
            pause(0.01)
            if strcmp(get(gco,'Type'),'uicontrol')
                if strcmp(get(gco,'String'),'done')

                    %-Get d allowing any 2nd order effects
                    %----------------------------------------------------------
                    for i = 1:m
                        for j = 1:m
                            if a(i,j)==1
                                d(i,j,k) = get(h4(i,j),'Value');
                            end
                        end
                    end
                    delete([h4(find(a==1))])
                    break

                end
            end
        end
    end
    delete(q)
else
    nlDCM = 0;
end
            
            
%==========================================================================
% slice timing
%==========================================================================
delays = spm_input('Slice timings [s]', -1, 'r', SPM.xY.RT*ones(1, m), m, [0 SPM.xY.RT]);


%==========================================================================
% echo time (TE) of data acquisition
%==========================================================================
TE    = 0;
TE_ok = 0;
while ~TE_ok
    TE = spm_input('Echo time, TE [s]');
    if ~TE || (TE < 0) || (TE > 0.1)
        str = { 'Extreme value for TE or TE undefined.',...
            'Please re-enter TE (note this value must be in seconds!).'};
        spm_input(str,1,'bd','OK',[1],1);
    else
        TE_ok = 1;
    end
end


spm_input('Thank you',1,'d')


%-Confounds (NB: the data have been filtered and whitened)
%--------------------------------------------------------------------------
v     = size(xY(1).u,1);
X0    = xY(1).X0;

%-Response variables
%--------------------------------------------------------------------------
n     = length(xY);
Y.dt  = SPM.xY.RT;
Y.X0  = X0;
for i = 1:n
    % regional responses
    Y.y(:,i)  = xY(i).u;
    Y.name{i} = xY(i).name;
end

%-Error precision components (one for each region) - i.i.d. (because of W)
%--------------------------------------------------------------------------
Y.Q        = spm_Ce(ones(1,n)*v);

%-Store all variables in DCM structure
%--------------------------------------------------------------------------
DCM.a      = a;
DCM.b      = b;
DCM.c      = c;
if nlDCM
    DCM.d  = d;
end
DCM.U      = U;
DCM.Y      = Y;
DCM.xY     = xY;
DCM.v      = v;
DCM.n      = n;
DCM.delays = delays;
DCM.TE     = TE;

%-Save
%--------------------------------------------------------------------------
if spm_matlab_version_chk('7') >= 0
    save(fullfile(swd,['DCM_' name]),'-V6','DCM');
else
    save(fullfile(swd,['DCM_' name]),'DCM');
end;
return
