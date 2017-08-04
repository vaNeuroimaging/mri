function Batch_wb_image_capture(ciftifile,outname,image_view,cifti_inds,palette,disp_zero,colorscaling,colorthresholding,Lsurf,Rsurf,Lsulc,Rsulc)
%Batch_wb_image_capture(ciftifile,outname,[image_view],[cifti_inds],[palette],[disp_zeros],[colorscaling],[colorthresholding],[Lsurf],[Rsurf],[Lsulc],[Rsulc])
%
% Automatically capture images of cifti data on surfaces without having to
% use the Workbench GUI. Requires wb_command in your linux path and the
% /home/data/scripts/Resources/ folder (and subfolders)  in your matlab path.
%
% Inputs:
%
%  ciftifile - data to capture. Should be a multi-column dtseries.nii or
%   dscalar.nii file
%
%  outname - root name of output image(s). Output files will be named
%   [outname]_map[#]_[image_view].png
%
%  image_view - hemisphere (Left, Right) and view (Lateral, Medial,
%   Anterior, Posterior, Dorsal, Ventral) to be shown in captured images.
%   Leave blank ([]) or omit for default view, or specify one of: 
%   'LR_LM' (default), 'LR_L','LR_M','L_LM','R_LM','LR_A','LR_P','LR_D','LR_V','L_L','L_R','R_L','R_M','L_A','L_P','L_D','L_V','R_A','R_P','R_D','R_V'
%
%  cifti_inds - a vector specifying which columns of the input cifti to
%   capture images of. Leave blank ([]) or omit to capture images of all columns.
%
%  palette - color palette to use when mapping values in cifti. Leave blank ([])
%   or omit to use the default palette (or the palette already saved in the
%   cifti file, if any), or specify one of: 
%   'PSYCH','PSYCH-NO-NONE','ROY-BIG','ROY-BIG-BL','Orange-Yellow','Gray_Interp_Positive','Gray_Interp','clear_brain','videen_style','fidl','raich4_clrmid','raich6_clrmid','HSB8_clrmid','RBGYR20','RBGYR20P','POS_NEG','red-yellow','blue-lightblue','FSL','power_surf','fsl_red','fsl_green','fsl_blue','fsl_yellow','JET256'
%
%  disp_zero - specify whether to display zero values (set to 1) or not
%   (set to 0, omit, or leave blank ([])). 
%
%  colorscaling - a 1x4 vector specifying the scaling of the color palette.
%   The vector should be formatted as [POS_MIN POS_MAX NEG_MIN NEG_MAX],
%   where entries represent numeric values in the cifti data. Leave blank ([]) or
%   omit to use auto-scaling.
%
%  colorthresholding - a 1x4 vector specifing how to threshold values in
%   the cifti file. The vector should be formatted as 
%   [MIN_INSIDE_THRESHOLD  MAX_INSIDE_THRESHOLD  MIN_OUTSIDE_THRESHOLD  MAX_OUTSIDE_THRESHOLD].
%   Either the INSIDE or OUTSIDE values should be set to NaN, indicating to
%   use the type of thresholding with real (non-NaN) values. With INSIDE
%   thresholding, only values between MIN_INSIDE_THRESHOLD and
%   MAX_INSIDE_THRESHOLD will be displayed; with OUTSIDE thresholding, only
%   values lower than MIN_OUTSIDE_THRESHOLD and higher than
%   MAX_OUTSIDE_THRESHOLD will be displayed. Leave blank or omit to display
%   all values.
%
%  Lsurf - path to the left hemisphere surface to use for display. Leave blank ([]) or
%   omit to use the default fs_LR_32k inflated surface.
%
%  Rsurf - path to the right hemisphere surface to use for display. Leave blank ([]) or
%   omit to use the default fs_LR_32k inflated surface.
%
%  Lsulc - path to the left hemisphere sulcal depth map to display as an underlay
%   beneath the cifti data. Specify 'none' to use no sulcal depth underlay.
%   Leave blank ([]) or omit to use the default Conte69 average sulcal depth map.
%
%  Rsulc - path to the right hemisphere sulcal depth map to display as an underlay
%   beneath the cifti data. Specify 'none' to use no sulcal depth underlay.
%   Leave blank ([]) or omit to use the default Conte69 average sulcal depth map.
%
%
%E.Gordon 06/17

origcapturefolder = '/home/data/scripts/Resources/batch_image_capture/';

[path,~,~] = fileparts(outname);
if isempty(path)
    path = pwd;
end

capturefolder = [path '/temp_image_capture_files/'];
mkdir(capturefolder)
copyfile([origcapturefolder '/*'],capturefolder)


if ~exist('image_view') || isempty(image_view)
    image_view = 'LR_LM';
end
if ~exist('Lsurf') || isempty(Lsurf)
    Lsurf = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.inflated.32k_fs_LR.surf.gii';
end
copyfile(Lsurf,[capturefolder '/L.surf.gii']);

if ~exist('Rsurf') || isempty(Rsurf)
    Rsurf = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.inflated.32k_fs_LR.surf.gii';
end
copyfile(Rsurf,[capturefolder '/R.surf.gii']);

if ~exist('Lsulc') || isempty(Lsulc)
    Lsulc = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sulc.32k_fs_LR.shape.gii';
end
if strcmp(Lsulc,'none')
    temp = gifti([capturefolder '/L.surf.gii']);
    save(gifti(single(zeros(size(temp.vertices,1),1))),[capturefolder '/L.sulc.shape.gii'])
else
    copyfile(Lsulc,[capturefolder '/L.sulc.shape.gii']);
end

if ~exist('Rsulc') || isempty(Rsulc)
    Rsulc = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sulc.32k_fs_LR.shape.gii';
end
if strcmp(Lsulc,'none')
    temp = gifti([capturefolder '/R.surf.gii']);
    save(gifti(single(zeros(size(temp.vertices,1),1))),[capturefolder '/R.sulc.shape.gii'])
else
    copyfile(Lsulc,[capturefolder '/R.sulc.shape.gii']);
end


data = ft_read_cifti_mod(ciftifile);
if ~exist('cifti_inds') || isempty(cifti_inds)
    cifti_inds = [1:size(data.data,2)];
end

out = data;
for ind = cifti_inds(:)'
    out.data = data.data(:,ind);
    out.mapname = {'Column number'};
    out.dimord = 'scalar_pos';
    ft_write_cifti_mod([capturefolder 'cifti'],out)
    
    if (exist('palette') && ~isempty(palette)) || (exist('colorscaling') && ~isempty(colorscaling)) || (exist('colorthresholding') && ~isempty(colorthresholding)) || (exist('disp_zero') && ~isempty(disp_zero))
        wbstring = ['wb_command -cifti-palette ' capturefolder '/cifti.dscalar.nii'];
        
        if (exist('colorscaling') && ~isempty(colorscaling))
            wbstring = [wbstring ' MODE_USER_SCALE ' capturefolder '/cifti2.dscalar.nii -pos-user ' num2str(colorscaling(1)) ' ' num2str(colorscaling(2)) ' -neg-user ' num2str(colorscaling(3)) ' ' num2str(colorscaling(4))];
        else
            wbstring = [wbstring ' MODE_AUTO_SCALE_PERCENTAGE ' capturefolder '/cifti2.dscalar.nii'];
        end
        if (exist('palette') && ~isempty(palette))
            wbstring = [wbstring ' -palette-name ' palette];
        end
        if (exist('colorthresholding') && ~isempty(colorthresholding))
            if all(isnan(colorthresholding([3:4])))
                wbstring = [wbstring ' -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_INSIDE ' num2str(colorthresholding(1)) ' ' num2str(colorthresholding(2))];
            elseif all(isnan(colorthresholding([1:2])))
                wbstring = [wbstring ' -thresholding THRESHOLD_TYPE_NORMAL THRESHOLD_TEST_SHOW_OUTSIDE ' num2str(colorthresholding(3)) ' ' num2str(colorthresholding(4))];
            else
                error('either inside or outside thresholds must have NaN values')
            end
        end
        if (exist('disp_zero') && ~isempty(disp_zero))
            disp_zero_str = 'false';
            if disp_zero
               disp_zero_str = 'true';
            end
            wbstring = [wbstring ' -disp-zero ' disp_zero_str];
        end
        [~,~] = system(wbstring);
        copyfile([capturefolder '/cifti2.dscalar.nii'],[capturefolder '/cifti.dscalar.nii'])
        delete([capturefolder '/cifti2.dscalar.nii']);
    end
    
    height = 925;
    width = 1401;
    if any(strcmp(image_view,{'LR_L','LR_M','L_LM','R_LM'})) 
        width = ceil(width ./ 2);
    end
    
    [~,~] = system(['wb_command -show-scene ' capturefolder '/Capture.scene ' image_view ' ' outname '_map' num2str(ind) '_' image_view '.png ' num2str(width) ' ' num2str(height)]);
    
end

rmdir(capturefolder,'s')