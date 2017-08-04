function Batch_PPI_parametric

%USER INPUT (and more down below)
%--------------------------------------------------------------------------

%Names of subjects to be run.
subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274', 

%Full paths of masks to be used. A separate PPI will be run on each mask.
anatomicalmasks = {'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Andrew/Subjects/Connectivity/OneSampleT/Deluca_PCC_FirstRest/dACC_roi.mat'};

%Combinations of conditions (i.e. column numbers) within the SPM design matrix that will interact with eigenvariate. A separate PPI will be run on each condition combination.
designmatrixcolumns_toconvolvewith = {[1]};

%For each condition combination, a vector of the weights given to each condition.
%The order of vectors here is the same as the order of condition combinations in "designmatrixcolumns_toconvolvewith" above
designmatrixcolumn_weights = {[1]};

%Names of each condition combination
%The order of names here is the same as the order of condition combinations in "designmatrixcolumns_toconvolvewith" above
Conditionnames = {'Parametric'};

%TR of the functional data
TR = 2;

%Location of PPI Batch Template file
PPItemplatefile = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/PPI_Batch_Template.mat';

%The prefix of any additional preprocessing steps that were done to the data to be used in the PPI, that were not done to the data in the SPM analysis. 
%If none, provide an empty matrix []. 
additionalpreprocessingprefix_inPPI = [];

%Defining the ROI:
% 1 = The ROI is defined as a sphere around the peak responsive voxel within the anatomical mask, clipped by that mask.
% 0 = The ROI is defined as the anatomical mask
usecontrastpeakswithinmask = 0;


%The following 3 variables are only relevant if "usecontrastpeakswithinmask = 1" above. Otherwise, use contrast_number=1, pval_forcontrast=1, and radius=anything 

%Contrast used to define the peak voxel
contrast_number = 1;
%P-value used to clip active voxels
pval_forcontrast = 1;
%Radius of sphere around peak voxel
radius = 4;

%END USER INPUT
%--------------------------------------------------------------------------

conditions = {'1Back','2Back','3Back'};
condtimes = {[48 92 136] [4 114 158] [26 70 180]};
conddur = 12;
nbacktimepointstouse = [];
for cond = 1:length(conditions)
    condtimepointstouse{cond} = [condtimes{cond}(1):condtimes{cond}(1)+conddur-1 , condtimes{cond}(2):condtimes{cond}(2)+conddur-1, condtimes{cond}(3):condtimes{cond}(3)+conddur-1];
    nbacktimepointstouse = [nbacktimepointstouse condtimepointstouse{cond}];
end
nbacktimepointstouse = sort(nbacktimepointstouse);


currentdirectory = pwd;

%Subject loop
for subjectnum = 1:length(subjects)
    
    %name of this subject
    subjid = subjects{subjectnum};
    
    
    %USER INPUT (and more down below)
    %--------------------------------------------------------------------------
    
    %Location of SPM analysis
    SPMfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/SPM8/Parametric/'];
    
    %END USER INPUT
    %--------------------------------------------------------------------------
    
    SPMfilename = [SPMfolder 'SPM.mat'];
    
    %get needed info for this SPM analysis and for the defined contrast
    [SPM,xSPM] = spm_getSPM_PPI(SPMfilename,contrast_number,pval_forcontrast);

    %Loop through masks
    for masknum = 1:length(anatomicalmasks)
        
        %name of this mask
        anatomicalmaskname = anatomicalmasks{masknum};
        
        %get the "short" name for the current mask, i.e. without the path and file extension
        slashlocations = findstr(anatomicalmaskname,'/'); roilocation = findstr(anatomicalmaskname,'_roi.mat');
        shortmaskname = anatomicalmaskname(slashlocations(end)+1:roilocation(end)-1);
        
        %Loop through conditions
        for condition = 1:length(designmatrixcolumns_toconvolvewith)
            
            %Loop through sessions
            for session = 1:length(SPM.Sess)
                
                %Figure out the correct column in the design matrix for the current condition
                designmatrixcolumns = designmatrixcolumns_toconvolvewith{condition}';
                
                %Define the name of the PPI as "PPI_maskname_conditionname_sessionnumber"
                if length(SPM.Sess) > 1
                    PPIname = [shortmaskname '_' Conditionnames{condition} '_Sess' num2str(session)];
                else
                    PPIname = [shortmaskname '_' Conditionnames{condition}];
                end
                
                
                %USER INPUT
                %--------------------------------------------------------------------------
                
                %Define the PPI output folder
                PPIoutputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/PPI/' PPIname '/'];
                
                %END USER INPUT
                %--------------------------------------------------------------------------
                
                
                %Make the PPI output folder, overwriting any previous results
                try;rmdir(PPIoutputfolder,'s');catch;end
                mkdir(PPIoutputfolder);
                
                %If you want to make an ROI around the peak voxel
                if usecontrastpeakswithinmask
                    
                    roi = maroi('load', anatomicalmaskname);
                    sp = mars_space([SPMfolder 'con_0001.hdr']);
                    save_as_image(roi, [PPIoutputfolder shortmaskname '.nii'], sp);
                    
                    
                    %Do a bunch of stuff to figure out which voxels are within the anatomical mask, and find the peak within those voxels
                    Anat.spec = spm_vol([PPIoutputfolder shortmaskname '.nii']);
                    Q = ones(1,size(xSPM.XYZmm,2));
                    XYZ = Anat.spec.mat \ [xSPM.XYZmm; Q];
                    j  = find(spm_sample_vol(Anat.spec, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
                    xSPMwithinmaskindices = zeros(size(xSPM.Z));
                    xSPMwithinmaskindices(j) = 1;
                    xSPMZstemp = xSPM.Z .* xSPMwithinmaskindices;
                    
                    %Get the peak value
                    maxZval = max(xSPMZstemp(find(xSPMZstemp~=0)));
                    
                    %Use that to find the peak voxel coordinates
                    maxcoords = xSPM.XYZmm(:,find(xSPMZstemp==maxZval));
                    
                    clear xSPMZstemp maxZval Q XYZ j xSPMwithinmaskindices Anat
                    
                    
                    %Make the ROI sphere and write it out as an roi
                    roi_params.centre = maxcoords';
                    roi_params.radius = radius;
                    [o, others] = maroi_sphere(roi_params);
                    o = descrip(o,'sphere');
                    o = label(o,'sphere');
                    saveroi(o, [PPIoutputfolder 'sphere_roi.mat']);
                    
                    %Trim the ROI sphere with the anatomical mask and write it out as an roi
                    rois_to_combine = [sprintf('%-300s',[PPIoutputfolder 'sphere_roi.mat']);sprintf('%-300s',anatomicalmaskname)];
                    roilist = maroi('load_cell', rois_to_combine);
                    func = 'r1 & r2';
                    for i = 1:length(roilist)
                        eval(sprintf('r%d = roilist{%d};', i, i));
                    end
                    eval(['o=' func ';']);
                    o = label(o, func);
                    saveroi(o, [PPIoutputfolder 'trimmedsphere_roi.mat']);
                    
                    %Export the trimmed sphereical roi as an image in the same space as the con_0001 image of the SPM analysis 
                    roi = maroi('load', [PPIoutputfolder 'trimmedsphere_roi.mat']);
                    sp = mars_space([SPMfolder 'con_0001.hdr']);
                    save_as_image(roi, [PPIoutputfolder 'trimmedsphere.nii'], sp);
                    
                    %Define the mask to use as that trimmed sphere image
                    maskname = [PPIoutputfolder 'trimmedsphere.nii'];
                    
                else
                    
                    %Otherwise, export the anatomical mask as an image in the same space as the con_0001 image of the SPM analysis 
                    roi = maroi('load', anatomicalmaskname);
                    sp = mars_space([SPMfolder 'con_0001.hdr']);
                    save_as_image(roi, [PPIoutputfolder shortmaskname '.nii'], sp);
                    
                    %And define the mask to use as the given mask
                    maskname = [PPIoutputfolder shortmaskname '.nii'];
                    maxcoords = [0;0;0];
                    
                end
                
                
                %Set up important stuff in the xY variable, which defines the mask used to get the eigenvariate
                xY.name = PPIname;
                xY.Ic = 0;
                xY.Sess = session;
                
                %Get the eigenvariate
                [Y,xY] = spm_regions_PPI(xSPM, SPM, maxcoords, xY, maskname);
                
                %Define the conditions and weightings used in the PPI
                weightingmatrix = [designmatrixcolumns ones(length(designmatrixcolumns),1)*2 designmatrixcolumn_weights{condition}'];
                
                %Convolve the eigenvariate with the psychological conditions
                PPI = spm_peb_ppi_parametric(SPMfilename,'ppi',xY,weightingmatrix,PPIname,1);
                
                %Load the PPI template
                load(PPItemplatefile);
                
                %Put the TR into the template
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
                
                %Put the ouput folder into the template
                matlabbatch{1}.spm.stats.fmri_spec.dir{1} = PPIoutputfolder;
                
                %Figure out which scans from the original SPM analysis were in the current session, and load those scans into the template
                
                for scan = 1:length(nbacktimepointstouse)
                    fullscanname = SPM.xY.P(SPM.Sess(session).row(nbacktimepointstouse(scan)),:);
                    slashlocations = find(fullscanname=='/');
                    newscanname = [fullscanname(1:slashlocations(end)) additionalpreprocessingprefix_inPPI fullscanname(slashlocations(end)+1:end)];
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{scan} = strtrim(newscanname);
                end
                
                %Load the three terms (psychological condition, eigenvariate, and interaction of the two) into the template
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI-interaction';
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val  = PPI.ppi(nbacktimepointstouse);
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = Conditionnames{condition};
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val  = PPI.P(nbacktimepointstouse);
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'Eigenvariate';
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val  = PPI.Y(nbacktimepointstouse);
                
                %Save and run the template
                save([PPIoutputfolder 'PPIBatch.mat'], 'matlabbatch');
                spm_jobman('run',[PPIoutputfolder 'PPIBatch.mat']);
                
                clear matlabbatch
                
                
                
            end
            
        end
    end
    cd(currentdirectory);
end



end



function [SPM,xSPM] = spm_getSPM_PPI(spmmatfile,contrastnum,pvalthresh)


%-Select SPM.mat & note SPM results directory
%-----------------------------------------------------------------------

Ic = contrastnum;
pvalue = pvalthresh;

sts = 1;

%-Preliminaries...
%=======================================================================

%-Load SPM.mat
%-----------------------------------------------------------------------
swd = spm_str_manip(spmmatfile,'H');
try
    load(fullfile(swd,'SPM.mat'));
catch
    error(['Cannot read ' fullfile(swd,'SPM.mat')]);
end
SPM.swd = swd;


%-Change directory so that relative filenames are valid
%----------------------------------------------------------------------
cd(SPM.swd);

%-----------------------------------------------------------------------

xX   = SPM.xX;				    %-Design definition structure
XYZ  = SPM.xVol.XYZ;			%-XYZ coordinates
S    = SPM.xVol.S;			    %-search Volume {voxels}
R    = SPM.xVol.R;			    %-search Volume {resels}
M    = SPM.xVol.M(1:3,1:3);		%-voxels to mm matrix
VOX  = sqrt(diag(M'*M))';		%-voxel dimensions
titlestr = SPM.xCon(Ic).name;



%-Contrast definitions
%=======================================================================

%-Load contrast definitions (if available)
%-----------------------------------------------------------------------
try
    xCon = SPM.xCon;
catch
    xCon = {};
end


%=======================================================================
% - C O N T R A S T S ,  S P M   C O M P U T A T I O N ,   M A S K I N G
%=======================================================================

%-Get contrasts
%-----------------------------------------------------------------------


xCon = SPM.xCon;

IcAdd = [];

nc       = length(Ic);  % Number of contrasts

if (nc > 1)
    if nc==2
        But='Conjunction|Global';      Val=[1 nc];
    else
        But='Conj''n|Intermed|Global'; Val=[1 NaN nc];
    end
    n = spm_input('Null hyp. to assess?','+1','b',But,Val,1);
    if isnan(n)
        if nc==3,
            n = nc-1;
        else
            n = nc-spm_input('Effects under null ','0','n1','1',nc-1);
        end
    end
else
    n = 1;
end

%-Get contrasts for masking
%-----------------------------------------------------------------------

maskCon = 0;

Im = [];
pm = [];
Ex = [];


%-Create/Get title string for comparison
%-----------------------------------------------------------------------
str  = xCon(Ic).name;
mstr = 'masked [incl.] by';




%-Compute & store contrast parameters, contrast/ESS images, & SPM images
%=======================================================================
SPM.xCon = xCon;
if ~isfield(SPM.xX, 'fullrank')
    SPM  = spm_contrasts(SPM, unique([Ic, Im, IcAdd]));
else
    SPM  = spm_eeg_contrasts_conv(SPM, unique([Ic, Im, IcAdd]));
    xSPM = [];
    return;
end

xCon     = SPM.xCon;
STAT     = xCon(Ic(1)).STAT;
VspmSv   = cat(1,xCon(Ic).Vspm);



%-Degrees of Freedom and STAT string describing marginal distribution
%-----------------------------------------------------------------------
df          = [xCon(Ic(1)).eidf xX.erdf];
str = '';

switch STAT
    case 'T'
        STATstr = sprintf('%c%s_{%.0f}','T',str,df(2));
    case 'F'
        STATstr = sprintf('%c%s_{%.0f,%.0f}','F',str,df(1),df(2));
    case 'P'
        STATstr = sprintf('%s^{%0.2f}','PPM',df(1));
end



%-Compute conjunction as minimum of SPMs
%-----------------------------------------------------------------------
Z         = Inf;
for i     = Ic
    Z = min(Z,spm_get_data(xCon(i).Vspm,XYZ));
end

% P values for False Discovery FDR rate computation (all search voxels)
%=======================================================================
switch STAT
    case 'T'
        Ps = (1 - spm_Tcdf(Z,df(2))).^n;
    case 'P'
        Ps = (1 - Z).^n;
    case 'F'
        Ps = (1 - spm_Fcdf(Z,df)).^n;
end


%-Compute mask and eliminate masked voxels
%-----------------------------------------------------------------------
for i = Im
    
    Mask = spm_get_data(xCon(i).Vspm,XYZ);
    um   = spm_u(pm,[xCon(i).eidf,xX.erdf],xCon(i).STAT);
    if Ex
        Q = Mask <= um;
    else
        Q = Mask >  um;
    end
    XYZ       = XYZ(:,Q);
    Z         = Z(Q);
    if isempty(Q)
        fprintf('\n')                                           %-#
        warning(sprintf('No voxels survive masking at p=%4.2f',pm))
        break
    end
end

%-clean up interface
%-----------------------------------------------------------------------


%=======================================================================
% - H E I G H T   &   E X T E N T   T H R E S H O L D S
%=======================================================================

%-Height threshold - classical inference
%-----------------------------------------------------------------------
u      = -Inf;
k      = 0;
if STAT ~= 'P'
    
    
    %-Get height threshold
    %-------------------------------------------------------------------
    thresDesc = 'none';
    
    switch thresDesc
        
        
        
        case 'FWE' % family-wise false positive rate
            %---------------------------------------------------------------
            try
                u=xSPM.u;
            catch
                u  = spm_input('p value (family-wise error)','+0','r',0.05,1,[0,1]);
            end;
            thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
            u  = spm_uc(u,df,STAT,R,n,S);
            
        case 'FDR' % False discovery rate
            %---------------------------------------------------------------
            try
                u=xSPM.u;
            catch
                u  = spm_input('p value (false discovery rate)','+0','r',0.05,1,[0,1]);
            end;
            thresDesc = ['p<' num2str(u) ' (' thresDesc ')'];
            u  = spm_uc_FDR(u,df,STAT,n,VspmSv,0);
            
        otherwise  %-NB: no adjustment
            % p for conjunctions is p of the conjunction SPM
            %---------------------------------------------------------------
            u = pvalue;
            if u <= 1
                thresDesc = ['p<' num2str(u) ' (unc.)'];
                u = spm_u(u^(1/n),df,STAT);
            else
                thresDesc = [STAT '=' num2str(u) ];
            end
            
    end
    
    %-Height threshold - Bayesian inference
    %-----------------------------------------------------------------------
end % (if STAT)

%-Calculate height threshold filtering
%-------------------------------------------------------------------
Q      = find(Z > u);

%-Apply height threshold
%-------------------------------------------------------------------
Z      = Z(:,Q);
XYZ    = XYZ(:,Q);
if isempty(Q)
    warning(sprintf('No voxels survive height threshold u=%0.2g',u))
end


%-Extent threshold (disallowed for conjunctions)
%-----------------------------------------------------------------------
if ~isempty(XYZ) & nc == 1
    
    %-Get extent threshold [default = 0]
    %-------------------------------------------------------------------
    k = 0;
    %-Calculate extent threshold filtering
    %-------------------------------------------------------------------
    A     = spm_clusters(XYZ);
    Q     = [];
    for i = 1:max(A)
        j = find(A == i);
        if length(j) >= k; Q = [Q j]; end
    end
    
    % ...eliminate voxels
    %-------------------------------------------------------------------
    Z     = Z(:,Q);
    XYZ   = XYZ(:,Q);
    if isempty(Q)
        warning(sprintf('No voxels survive extent threshold k=%0.2g',k))
    end
    
else
    
    k = 0;
    
end % (if ~isempty(XYZ))


%=======================================================================
% - E N D
%=======================================================================


%-Assemble output structures of unfiltered data
%=======================================================================
xSPM   = struct( ...
    'swd',		swd,...
    'title',	titlestr,...
    'Z',		Z,...
    'n',		n,...
    'STAT',		STAT,...
    'df',		df,...
    'STATstr',	STATstr,...
    'Ic',		Ic,...
    'Im',		Im,...
    'pm',		pm,...
    'Ex',		Ex,...
    'u',		u,...
    'k',		k,...
    'XYZ',		XYZ,...
    'XYZmm',	SPM.xVol.M(1:3,:)*[XYZ; ones(1,size(XYZ,2))],...
    'S',		SPM.xVol.S,...
    'R',		SPM.xVol.R,...
    'FWHM',		SPM.xVol.FWHM,...
    'M',		SPM.xVol.M,...
    'iM',		SPM.xVol.iM,...
    'DIM',		SPM.xVol.DIM,...
    'VOX',		VOX,...
    'Vspm',		VspmSv,...
    'Ps',		Ps,...
    'thresDesc',thresDesc);

% RESELS per voxel (density) if it exists
%-----------------------------------------------------------------------
try, xSPM.VRpv  = SPM.VRpv;       end
try, xSPM.units = SPM.xVol.units; end

end




function [Y,xY] = spm_regions_PPI(xSPM,SPM,xyz,xY,maskname)

noGraph = 1;


xY.xyz = xyz;


%-Specify VOI
%--------------------------------------------------------------------------
xY.M = xSPM.M;
xY.spec = spm_vol(maskname);
xY.def = 'mask';
[xY, xY.XYZmm, Q] = spm_ROI(xY, xSPM.XYZmm);
% Q = ones(1,size(xSPM.XYZmm,2));
% XYZ    = xY.spec.mat \ [xSPM.XYZmm; Q];
% j  = find(spm_sample_vol(xY.spec, XYZ(1,:), XYZ(2,:), XYZ(3,:),0) > 0);
% xY.XYZmm = xSPM.XYZmm(:,j);
try, xY = rmfield(xY,'M'); end
try, xY = rmfield(xY,'rej'); end

if isempty(xY.XYZmm), error('Empty region.'); end

%-Extract required data from results files
%==========================================================================
spm('Pointer','Watch')

%-Get raw data, whiten and filter
%--------------------------------------------------------------------------
y        = spm_get_data(SPM.xY.VY,xSPM.XYZ(:,Q));
y        = spm_filter(SPM.xX.K,SPM.xX.W*y);


%-Computation
%==========================================================================

%-Remove null space of contrast
%--------------------------------------------------------------------------
if xY.Ic
    
    %-Parameter estimates: beta = xX.pKX*xX.K*y
    %----------------------------------------------------------------------
    beta  = spm_get_data(SPM.Vbeta,xSPM.XYZ(:,Q));
    
    %-subtract Y0 = XO*beta,  Y = Yc + Y0 + e
    %----------------------------------------------------------------------
    y     = y - spm_FcUtil('Y0',SPM.xCon(xY.Ic),SPM.xX.xKXs,beta);
    
end

%-Confounds
%--------------------------------------------------------------------------
xY.X0     = SPM.xX.xKXs.X(:,[SPM.xX.iB SPM.xX.iG]);

%-Extract session-specific rows from data and confounds
%--------------------------------------------------------------------------
try
    i     = SPM.Sess(xY.Sess).row;
    y     = y(i,:);
    xY.X0 = xY.X0(i,:);
end

% and add session-specific filter confounds
%--------------------------------------------------------------------------
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).X0];
end
try
    xY.X0 = [xY.X0 SPM.xX.K(xY.Sess).KH]; % Compatibility check
end

%-Remove null space of X0
%--------------------------------------------------------------------------
xY.X0     = xY.X0(:,any(xY.X0));


%-Compute regional response in terms of first eigenvariate
%--------------------------------------------------------------------------
[m n]   = size(y);
if m > n
    [v s v] = svd(y'*y);
    s       = diag(s);
    v       = v(:,1);
    u       = y*v/sqrt(s(1));
else
    [u s u] = svd(y*y');
    s       = diag(s);
    u       = u(:,1);
    v       = y'*u/sqrt(s(1));
end
d       = sign(sum(v));
u       = u*d;
v       = v*d;
Y       = u*sqrt(s(1)/n);

%-Set in structure
%--------------------------------------------------------------------------
xY.y    = y;
xY.u    = Y;
xY.v    = v;
xY.s    = s;
end




function PPI = spm_peb_ppi_parametric(varargin)
% Bold deconvolution to create physio- or psycho-physiologic interactions
% FORMAT PPI = spm_peb_ppi(SPMname,ppiflag,VOI,Uu,ppiname,showGraphics)
%
% SPM          - Structure containing generic details about the analysis or
%                the fully qualified filename of such a structure.
% ppiflag      - Type of analysis. Must be one of:
%                  'simple deconvolution'          or 'sd'
%                  'psychophysiologic interaction' or 'ppi'
%                  'physiophysiologic interaction' or 'phipi'
% VOI          - Structure containing details about a VOI (as produced by
%                spm_regions) or the fully qualified filename of such a
%                structure. If a structure, then VOI should be of size 1x1
%                in the case of simple deconvolution, and psychophysiologic 
%                interactions) or 1x2, in the case of physiophysiologic
%                interactions. If a file name it should be 1xN or 2xN.
% Uu           - Matrix of input variables and contrast weights. This is an
%                [n x 3] matrix. The first column indexes SPM.Sess.U(i). The
%                second column indexes the name of the input or cause, see
%                SPM.Sess.U(i).name{j}. The third column is the contrast
%                weight. Unless there are parametric effects the second
%                column will generally be a 1.
% ppiname      - Basename of the PPI file to save. The saved file will be:
%                <PATH_TO_SPM.MAT>/PPI_<ppiname>.mat
% showGraphics - empty or 1 = yes, 0 = no.
%
%
% PPI.ppi    - (PSY*xn  or xn1*xn2) convolved with the HRF
% PPI.Y      - Original BOLD eigenvariate. Use as covariate of no interest
% PPI.P      - PSY convolved with HRF for psychophysiologic interactions,
%              or in the case of physiophysologic interactions contains
%              the eigenvariate of the second region. 
% PPI.name   - Name of PPI
% PPI.xY     - Original VOI information
% PPI.xn     - Deconvolved neural signal(s)
% PPI.U.u    - Psychological variable or input function (PPIs only)
% PPI.U.w    - Contrast weights for psychological variable (PPIs only)
% PPI.U.name - Names of psychological conditions (PPIs only)
%__________________________________________________________________________
%
% This routine is effectively a hemodynamic deconvolution using full priors
% and EM to deconvolve the HRF from a hemodynamic time series to give a 
% neuronal time series [that can be found in PPI.xn].  This deconvolution 
% conforms to Wiener filtering. The neuronal process is then used to form 
% PPIs. See help text within function for more details.
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Darren Gitelman
% $Id: spm_peb_ppi.m 4185 2011-02-01 18:46:18Z guillaume $

% SETTING UP A PPI THAT ACCOUNTS FOR THE HRF
% =========================================================================
% PPI's were initially conceived as a means of identifying regions whose
% reponses can be explained in terms of an interaction between activity in
% a specified source (the physiological factor) and some experimental
% effect (the psychological factor). However, a problem in setting up PPI's
% is that in order to derive a proper estimate of the interaction between
% a psychological variable (P) and measured hemodynamic signal (x), one 
% cannot simply convolve the psychological variable with the hrf (HRF) and 
% multiply by the signal. Thus:
% 
%                  conv(P,HRF).* x ~= conv((P.*xn),HRF)
%
% P   = psychological variable
% HRF = hemodynamic response function
% xn  = underlying neural signal which in fMRI is convolved with the hrf to
%       give the signal one measures -- x.
% x   = measured fmri signal
%
% It is actually the right hand side of the equation one wants.
% Thus one has to work backwards, in a sense, and deconvolve the hrf
% from x to get xn. This can then be multiplied by P and the resulting
% vector (or matrix) reconvolved with the hrf.
%
% This algorithm uses a least squares strategy to solve for xn.
%
% The source's hemodynamics are x = HRF*xn;
%
% Using the constraint that xn should have a uniform spectral density 
% we can expand x in terms of a discrete cosine set (xb)
%
%      xn  = xb*B
%       B  = parameter estimate
%
% The estimator of x is then
%
%       x  = HRF(k,:)*xn
%       x  = HRF(k,:) * xb * B
%
% This accounts for different time resolutions between our hemodynamic 
% signal and the discrete representation of the psychological variable. In 
% this case k is a vector representing the time resolution of the scans.
%
% Conditional estimates of B allow for priors that ensure uniform variance 
% over frequencies.
%
% PPI STATISTICAL MODEL
% =========================================================================
% Once the PPI.ppi interaction term has been calculated a new GLM must be
% setup to search for the interaction effects across the brain. This is
% done using a standard, first level, fMRI model, which must include 3
% covariates, PPI.ppi (interaction), PPI.Y (main effect: source region bold
% signal) and PPI.P (main effect: "psychological" condition), plus any
% nuisance regressors according to the particular design.
%
% NB: Designs that include only the interaction term without the main
% effects are not proper as inferences on the interaction will include a
% mixture of both main and interaction effects. 
%
% Once the model has been setup and run, a contrast of [1 0 0 ] over the
% PPI.ppi, PPI.Y and PPI.P columns respectively, will show regions with a
% positive relationship to the interaction term, discounting any main
% effects. Negative regressions can be examined with [-1 0 0]. A PPI random
% effects analysis would involve taking the con*.img files from the [1 0 0]
% t-contrast for each subject and forwarding them to a second level
% analysis.


% Set up the graphical interface
%--------------------------------------------------------------------------
Finter = spm_figure('GetWin','Interactive');
header = get(Finter,'Name');
spm_clf(Finter); set(Finter,'name','PPI Setup');


% Check inputs and set up variables
%--------------------------------------------------------------------------
if nargin>0 && isstruct(varargin{1})
    SPM = varargin{1};
else
    try
        load(varargin{1})
    catch
        [P, sts]     = spm_select(1,'^SPM\.mat$','Select SPM.mat');
        if ~sts, PPI = struct([]); return; end
        swd          = spm_str_manip(P,'H');
        load(fullfile(swd,'SPM.mat'));
        SPM.swd      = swd;
    end
end
cwd    = pwd;
cd(SPM.swd)

% Setup variables
%--------------------------------------------------------------------------
RT      = SPM.xY.RT;
dt      = SPM.xBF.dt;
NT      = RT/dt;
fMRI_T0 = SPM.xBF.T0;

% Ask whether to perform physiophysiologic or psychophysiologic interactions
%--------------------------------------------------------------------------
try
    ppiflag = varargin{2};
catch
    ppiflag = {'simple deconvolution',...
               'psychophysiologic interaction',...
               'physiophysiologic interaction'};
    i       = spm_input('Analysis type?',1,'m',ppiflag);
    ppiflag = ppiflag{i};
end


switch lower(ppiflag)
    
    case  {'simple deconvolution','sd'}
    %======================================================================
    if nargin>2 && isstruct(varargin{3})
        p.xY = varargin{3};
    else
        try
            VOI = varargin{3};
            p   = load(deblank(VOI(1,:)),'xY');
        catch
            spm_input('physiological variable:...  ',2,'d');
            voi = spm_select(1,'^VOI.*\.mat$',{'select VOI'});
            p   = load(deblank(voi),'xY');
        end
    end
    xY(1) = p.xY;
    Sess  = SPM.Sess(xY(1).Sess);

    case  {'physiophysiologic interaction','phipi'}
    %======================================================================
    if nargin>2 && isstruct(varargin{3})
        xY = varargin{3};
        xY = xY(:)';
        if size(xY) ~= [1 2]
            error('Must include 2 VOI structures for physiophysiologic interactions')
        end
    else
        try
            VOI = varargin{3};
            if size(VOI,1) ~= 2
                error('Must include 2 VOI filenames for physiophygiologic interactions')
            end
            for i = 1:2
                p     = load(deblank(VOI(i,:)),'xY');
                xY(i) = p.xY;
            end
        catch
            spm_input('physiological variables:...  ',2,'d');
            voi      = spm_select(2,'^VOI.*\.mat$',{'select VOIs'});
            for  i = 1:2
                p      = load(deblank(voi(i,:)),'xY');
                xY(i)  = p.xY;
            end
        end
    end
    Sess = SPM.Sess(xY(1).Sess);

    case  {'psychophysiologic interaction','ppi'}
    %======================================================================
    if nargin>2 && isstruct(varargin{3})
        p.xY = varargin{3};
    else
        try
            VOI = varargin{3};
            p   = load(deblank(VOI(1,:)),'xY');
        catch
            spm_input('physiological variable:...  ',2,'d');
            voi = spm_select(1,'^VOI.*\.mat$',{'select VOI'});
            p   = load(deblank(voi(:))','xY');
        end
    end
    xY(1) = p.xY;
    Sess  = SPM.Sess(xY(1).Sess);

    % get 'causes' or inputs U
    %----------------------------------------------------------------------
    U.name = {};
    U.u    = [];
    U.w    = [];
    try
        Uu = varargin{4};
        for i = 1:size(Uu,1)
            U.u           = [U.u Sess.U(Uu(i,1)).u(33:end,Uu(i,2))];
            U.name{end+1} = Sess.U(Uu(i,1)).name{Uu(i,2)};
            U.w           = [U.w Uu(i,3)];
        end
    catch
        spm_input('Psychological variable:...  ',2,'d');
        u      = length(Sess.U);
        for  i = 1:u
            for  j = 1:length(Sess.U(i).name)
                str   = ['include ' Sess.U(i).name{j} '?'];
                if spm_input(str,3,'y/n',[1 0])
                    str             = 'Contrast weight';
                    tmpw            = spm_input(str,4,'e',[],1);
                    % if tmpw==0 then don't include the column in the
                    % design. This takes care of the possibility that
                    % the user would select to include the column but
                    % then give it a 0 weight.
                    %------------------------------------------------
                    if tmpw ~= 0
                        U.w             = [U.w tmpw];
                        U.u             = [U.u Sess.U(i).u(33:end,j)];
                        U.name{end + 1} = Sess.U(i).name{j};
                    end
                end
            end
        end
    end

end % (switch setup)


% Name of PPI file to be saved
%--------------------------------------------------------------------------
try
    PPI.name = varargin{5};
catch
    PPI.name = spm_input('Name of PPI',3,'s','PPI');
end
[tmp ppiFilename] = fileparts(PPI.name);


% Check if Graphical output should be shown
%--------------------------------------------------------------------------
try
    showGraphics = varargin{6};
catch
    showGraphics = 1;
end
if showGraphics
    Fgraph       = spm_figure('GetWin','PPI');
    spm_clf(Fgraph);
    FS           = spm('FontSizes');
end


% Setup more variables
%--------------------------------------------------------------------------
N = length(xY(1).u);
k = 1:NT:N*NT;                             % microtime to scan time indices


% Create basis functions and hrf in scan time and microtime
%--------------------------------------------------------------------------
spm('Pointer','watch')
hrf = spm_hrf(dt);


% Create convolved explanatory {Hxb} variables in scan time
%--------------------------------------------------------------------------
xb  = spm_dctmtx(N*NT + 128,N);
Hxb = zeros(N,N);
for i = 1:N
    Hx       = conv(xb(:,i),hrf);
    Hxb(:,i) = Hx(k + 128);
end
xb = xb(129:end,:);


% Get confounds (in scan time) and constant term
%--------------------------------------------------------------------------
X0 = xY(1).X0;
M  = size(X0,2);


% Get response variable,
%--------------------------------------------------------------------------
for i = 1:size(xY,2)
    Y(:,i) = xY(i).u;
end


% Remove confounds and save Y in ouput structure
%--------------------------------------------------------------------------
Yc    = Y - X0*inv(X0'*X0)*X0'*Y;
PPI.Y = Yc(:,1);
if size(Y,2) == 2
    PPI.P = Yc(:,2);
end


% Specify covariance components; assume neuronal response is white
% treating confounds as fixed effects
%--------------------------------------------------------------------------
Q = speye(N,N)*N/trace(Hxb'*Hxb);
Q = blkdiag(Q, speye(M,M)*1e6  );


% Get whitening matrix (NB: confounds have already been whitened)
%--------------------------------------------------------------------------
W = SPM.xX.W(Sess.row,Sess.row);


% Create structure for spm_PEB
%--------------------------------------------------------------------------
clear P
P{1}.X = [W*Hxb X0];        % Design matrix for lowest level
P{1}.C = speye(N,N)/4;      % i.i.d assumptions
P{2}.X = sparse(N + M,1);   % Design matrix for parameters (0's)
P{2}.C = Q;


switch ppiflag

    case  {'simple deconvolution','sd'}
    %======================================================================
    C  = spm_PEB(Y,P);
    xn = xb*C{2}.E(1:N);
    xn = spm_detrend(xn);

    % Save variables (NOTE: xn is in microtime and does not account for
    % slice timing shifts). To convert to BOLD signal convolve with a hrf.
    % Use a microtime to scan time index to convert to scan time: e.g.,
    % k = 1:NT:N*NT; where NT = number of bins per TR = TR/dt or SPM.xBF.T
    % and N = number of scans in the session. Finally account for slice
    % timing effects by shifting the index accordingly. See approximately
    % lines 420 and 421 below for an example.)
    %----------------------------------------------------------------------
    PPI.xn = xn;

    % Plot results
    %----------------------------------------------------------------------
    if showGraphics
        figure(Fgraph);
        t = RT*[1:N];
        T = dt*[1:(N*NT)];

        ax = subplot(2,1,1);
        plot(t,Yc,T,PPI.xn)
        title('hemodynamic and neuronal responses')
        xlabel('time (secs)')
        axis tight square
        grid on
        legend('BOLD','neuronal')

        str = sprintf('Simple Deconvolution: %s\n',ppiFilename);
        str = [str sprintf('VOI file: %s',xY.name)];
        hAx = axes('Position',[0 0 1 1],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','Tex',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Parent',Fgraph,...
                'Units','points',...
                'Visible','off');
        AxPos = get(hAx,'Position'); 
        set(hAx,'XLim',[0,AxPos(3)]); set(hAx,'YLim',[0,AxPos(4)])
        dy    = FS(9);
        h     = text(dy,floor(AxPos(4))-2*dy,str);
        set(h,'Parent',hAx);
    end

    case  {'physiophysiologic interaction','phipi'}
    %======================================================================
    C    = spm_PEB(Y(:,1),P);
    xn1  = xb*C{2}.E(1:N);
    C    = spm_PEB(Y(:,2),P);
    xn2  = xb*C{2}.E(1:N);
    xn1  = spm_detrend(xn1);
    xn2  = spm_detrend(xn2);
    xnxn = xn1.*xn2;

    % Convolve, convert to scan time, and account for slice timing shift
    %----------------------------------------------------------------------
    ppi = conv(xnxn,hrf);
    ppi = ppi((k-1) + fMRI_T0);

    % Save variables
    %----------------------------------------------------------------------
    PPI.xn  = [xn1 xn2];
    PPI.ppi = spm_detrend(ppi);

    % Plot results
    %----------------------------------------------------------------------
    if showGraphics
        figure(Fgraph);
        t = RT*[1:N];
        T = dt*[1:(N*NT)];

        ax = subplot(2,1,1);
        plot(t,PPI.ppi)
        title('PPI')
        xlabel('time (secs)')
        axis tight square
        grid on

        subplot(2,2,3)
        plot(t,Yc(:,1),T,PPI.xn(:,1))
        title('hemodynamic and neuronal responses (1st)')
        xlabel('time (secs)')
        axis tight square
        grid on
        legend('BOLD','neuronal')

        subplot(2,2,4)
        plot(t,Yc(:,2),T,PPI.xn(:,2))
        title('hemodynamic and neuronal responses (2nd)')
        xlabel('time (secs)')
        axis tight square
        grid on
        legend('BOLD','neuronal')
        
        str = sprintf('Physiophysiologic Interaction: %s\n',ppiFilename);
        str = [str, sprintf('VOI File 1: %s\n',xY(1).name)];
        str = [str, sprintf('VOI File 2: %s',xY(2).name)];
        hAx = axes('Position',[0 0 1 1],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','Tex',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Parent',Fgraph,...
                'Units','points',...
                'Visible','off');
        AxPos = get(hAx,'Position'); 
        set(hAx,'XLim',[0,AxPos(3)]); set(hAx,'YLim',[0,AxPos(4)])
        dy    = FS(9);
        h     = text(dy,floor(AxPos(4))-2*dy,str);
        set(h,'Parent',hAx);
        
    end

    case  {'psychophysiologic interaction','ppi'}
    %======================================================================

    % COMPUTE PSYCHOPHYSIOLOGIC INTERACTIONS
    % use basis set in microtime
    %----------------------------------------------------------------------
    % get parameter estimates and neural signal; beta (C) is in scan time
    % This clever trick allows us to compute the betas in scan time which is
    % much quicker than with the large microtime vectors. Then the betas
    % are applied to a microtime basis set generating the correct neural
    % activity to convolve with the psychological variable in microtime
    %----------------------------------------------------------------------
    C  = spm_PEB(Y,P);
    xn = xb*C{2}.E(1:N);
    xn = spm_detrend(xn);

    % Setup psychological variable from inputs and contast weights
    %----------------------------------------------------------------------
    PSY = zeros(N*NT,1);
    for i = 1:size(U.u,2)
        PSY = PSY + full(U.u(:,i)*U.w(:,i));
    end
    % PSY = spm_detrend(PSY);  <- removed centering of psych variable
    % prior to multiplication with xn. Based on discussion with Karl
    % and Donald McLaren. 

    % Multiply psychological variable by neural signal
    %----------------------------------------------------------------------
    PSY = PSY+2;
    
    PSYxn = PSY.*xn;

    % Convolve, convert to scan time, and account for slice timing shift
    %----------------------------------------------------------------------
    ppi = conv(PSYxn,hrf);
    ppi = ppi((k-1) + fMRI_T0);

    % Convolve psych effect, convert to scan time, and account for slice
    % timing shift
    %----------------------------------------------------------------------
    PSYHRF = conv(PSY,hrf);
    PSYHRF = PSYHRF((k-1) + fMRI_T0);

    % Save psychological variables
    %----------------------------------------------------------------------
    PPI.psy = U;
    PPI.P   = PSYHRF;
    PPI.xn  = xn;
    PPI.ppi = spm_detrend(ppi);

    % Plot results
    %----------------------------------------------------------------------
    if showGraphics
        figure(Fgraph);
        t = RT*[1:N];
        T = dt*[1:(N*NT)];

        str = sprintf('Psychophyiologic Interaction: %s\n',ppiFilename);
        str = [str, sprintf('VOI File: %s\n',xY(1).name)];
        str = [str, sprintf('Factors: ')];
        for i = 1:numel(U.name)
            str = [str, sprintf('%s [%0.0f]',U.name{i},U.w(i))];
            if i < numel(U.name)
                str = [str, sprintf('; ')];
            end
        end

        ax = subplot(2,1,1);
        plot(t,Yc(:,1),T,PPI.xn(:,1))
        title('hemodynamic and neuronal responses')
        xlabel('time (secs)')
        axis tight square
        grid on
        legend('BOLD','neuronal')

        subplot(2,2,3)
        plot(T,PSY,'LineStyle','--','Color',[0 .65 0]);
        hold on
        plot(t,PPI.P,'LineStyle','-','LineWidth',1,'Color','b');
        hold off
        title('[convolved] psych. variable')
        xlabel('time (secs)')
        axis tight square
        grid on

        subplot(2,2,4)
        plot(t,PPI.ppi)
        title('PPI')
        xlabel('time (secs)')
        axis tight square
        grid on
        
        hAx = axes('Position',[0 0 1 1],...
                'DefaultTextFontSize',FS(8),...
                'DefaultTextInterpreter','None',...
                'DefaultTextVerticalAlignment','Baseline',...
                'Parent',Fgraph,...
                'Units','points',...
                'Visible','off');
        AxPos = get(hAx,'Position'); 
        set(hAx,'XLim',[0,AxPos(3)]); set(hAx,'YLim',[0,AxPos(4)])
        dy    = FS(9);
        h     = text(dy,floor(AxPos(4))-2*dy,str);
        set(h,'Parent',hAx);
    end
    
end % (switch)

% Setup other output variables and Save
%--------------------------------------------------------------------------
PPI.xY = xY;
PPI.dt = dt;
str    = ['PPI_' PPI.name '.mat'];

if spm_check_version('matlab','7') >= 0,
    save(fullfile(SPM.swd,str),'-V6','PPI')
else
    save(fullfile(SPM.swd,str),'PPI')
end

% Clean up
%--------------------------------------------------------------------------
spm('Pointer','arrow'); set(Finter,'name',header);
fprintf('\nCompleted PPI: %s\n',ppiFilename);
cd(cwd);

end





