function Batch_DCM

warning off

directory = pwd;

outputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DCM/DCMOutput.txt';

subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
    

regions = {'DC_bilat','raIns','Precuneus'};
models = {'Rest1','Rest2','Rest3','Rest4','Rest5','Rest6'};

regionfolder = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/DCM_ROIs/';

delete([outputfilename]);
fid = fopen([outputfilename],'at');
fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Subject','ROI_connection','Parameter','Model');
fclose(fid);
dlmwrite([outputfilename],' ','-append');

for modelnum = 1:length(models)
    model = models{modelnum};
    
    for subject = 1:length(subjects)
        
        disp(subjects{subject})
        
        datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/SPM8/FirstRest_forDCM/'];
        
        [SPM, xSPM] = spm_getSPM_Evan([datafolder 'SPM.mat']);
        
        for region = 1:length(regions)
            
            regionname = [regionfolder regions{region} '_roi.mat'];
            
            roi = maroi('load', regionname);
            sp = mars_space([datafolder 'con_0001.hdr']);
            save_as_image(roi, [datafolder regions{region} '.nii'], sp);
            
            %And set this saved image be the mask to be used
            imagename = [regions{region} '.nii'];
            maxcoords = [0;0;0];
            
            
            
            
            %Set up important stuff in the xY variable, which defines the mask used to get the eigenvariate
            xY.name = imagename;
            xY.Ic = 0;
            xY.Sess = 1;
            
            %Get the eigenvariate
            [Y,xY] = spm_regions_DCM(xSPM, SPM, maxcoords, xY, imagename);
            
            
            
        end
        
        
        load(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DCM/DCM_' model '.mat'])
        for region = 1:length(regions)
            load([datafolder 'VOI_' regions{region} '.mat']);
            
            DCM.Y.X0 = xY.X0;
            DCM.Y.y(:,region) = xY.u;
            DCM.Y.name{region} = regions{region};
            DCM.xY(region) = xY;
            
        end
        
        save([datafolder 'DCM_' model],'DCM');
        
        clear DCM
        
        spm_dcm_estimate([datafolder 'DCM_' model]);
        
        load([datafolder 'DCM_' model]);
        
        
        for i = 1:length(regions)
            for j = 1:length(regions)
                if i~=j
                    
                    texttowrite = [subjects{subject},'   ',[regions{i} '_to_' regions{j}],'   ',num2str(DCM.Ep.A(j,i)),'   ',model];
                    dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
                    
                end
            end
        end
        
        cd(directory)
         
        
%P(subject,:) = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/DCM_' model '.mat'];


    end
    
%     cd /fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DCM/
%     spm_dcm_average(1,P,['Group' model '.mat']);
%     clear P
     
    
end

end

function [Y,xY] = spm_regions_DCM(xSPM,SPM,xyz,xY,maskname)

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

%-Save
%==========================================================================
str = ['VOI_' xY.name(1:end-4) '.mat'];
% if isfield(xY,'Sess') && isfield(SPM,'Sess')
%     str = sprintf('VOI_%s_%i.mat',xY.name,xY.Sess);
% end
if spm_check_version('matlab','7') >= 0
    save(fullfile(SPM.swd,str),'-V6','Y','xY')
else
    save(fullfile(SPM.swd,str),'Y','xY')
end

fprintf('   VOI saved as %s\n',spm_str_manip(fullfile(SPM.swd,str),'k55'));

end




