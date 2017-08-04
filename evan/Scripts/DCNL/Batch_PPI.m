function Batch_PPI

%USER INPUT (and more down below)
%--------------------------------------------------------------------------

%Names of subjects to be run.
subjects = {'110','189','199','269'};
%{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
% 
%'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274', 

%Full paths of masks to be used. A separate PPI will be run on each mask.
anatomicalmasks = {'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/DC_bilat_roi.mat','/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/VSi_bilat_roi.mat'};

%Combinations of conditions (i.e. column numbers) within the SPM design matrix that will interact with eigenvariate. A separate PPI will be run on each condition combination.
designmatrixcolumns_toconvolvewith = {[1 2 3 4],[1 3]};

%For each condition combination, a vector of the weights given to each condition.
%The order of vectors here is the same as the order of condition combinations in "designmatrixcolumns_toconvolvewith" above
designmatrixcolumn_weights = {[1 1 1 -1],[-1 1]};

%Names of each condition combination
%The order of names here is the same as the order of condition combinations in "designmatrixcolumns_toconvolvewith" above
Conditionnames = {'Task_vs_Fix','3B_vs_1B'};

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

currentdirectory = pwd;

%Subject loop
for subjectnum = 1:length(subjects)
    
    %name of this subject
    subjid = subjects{subjectnum};
    
    
    %USER INPUT (and more down below)
    %--------------------------------------------------------------------------
    
    %Location of this subject's SPM analysis
    SPMfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/SPM8/Cond/'];
    
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
                    
                    %And set this saved image be the mask to be used
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
                weightingmatrix = [designmatrixcolumns ones(length(designmatrixcolumns),1) designmatrixcolumn_weights{condition}'];
                
                %Convolve the eigenvariate with the psychological conditions
                PPI = spm_peb_ppi(SPMfilename,'ppi',xY,weightingmatrix,PPIname,1);
                
                %Load the PPI template
                load(PPItemplatefile);
                
                %Put the TR into the template
                matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
                
                %Put the ouput folder into the template
                matlabbatch{1}.spm.stats.fmri_spec.dir{1} = PPIoutputfolder;
                
                %Figure out which scans from the original SPM analysis were in the current session, and load those scans into the template
                
                for scan = 1:length(SPM.Sess(session).row)
                    fullscanname = SPM.xY.P(SPM.Sess(session).row(scan),:);
                    slashlocations = find(fullscanname=='/');
                    newscanname = [fullscanname(1:slashlocations(end)) additionalpreprocessingprefix_inPPI fullscanname(slashlocations(end)+1:end)];
                    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{scan} = strtrim(newscanname);
                end
                
                %Load the three terms (psychological condition, eigenvariate, and interaction of the two) into the template
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).name = 'PPI-interaction';
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val  = PPI.ppi;
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).name = Conditionnames{condition};
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val  = PPI.P;
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).name = 'Eigenvariate';
                matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val  = PPI.Y;
                
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





