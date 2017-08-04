function map_vol_to_surface(volumepattern,hem,mappingtype,space)
%map_vol_to_surface(volume,[hem],[mappingtype],[space])
%
%     Uses workbench to map a nifti or 4dfp volume to a workbench surface so
%     that it can be displayed in workbench.
%
%     volume - Full path of nifti volume. The surface metric will be written
%     out in the same folder as the volume. Use of * indexing is supported.
%
%     hem - Hemispheres to map. Can be 'L', 'R', or 'both.' Defaults to
%     'both'.
%
%     mappingtype - Specify the method used to map the volume to the surface.
%      Can be 'enclosing','trilinear', or 'ribbon-constrained' (see below
%      for details). 'ribbon constrained' should generally be used for
%      scalar data; 'enclosing voxel' should be used for categorical data.
%      Defaults to 'ribbon-constrained'.
%
%    space - Space of input volume; can be 'MNI' or '711-2B'. Volumes must
%    be in MNI space to be mapped to the surface; if '711-2B' is selected,
%    the volume will be registered to MNI space before mapping. Defaults to
%    'MNI'. 
%
%
%    Mapping type descriptions (from the wb_command help text)
%         Enclosing voxel uses the value from the voxel the vertex lies
%         inside, while trilinear does a 3D linear interpolation based on
%         the voxels immediately on each side of the  vertex's position.
%         The ribbon mapping method constructs a polyhedron from the
%         vertex's neighbors on each surface, and estimates the amount of
%         this polyhedron's volume that falls inside any nearby voxels, to
%         use as the weights for sampling.  The volume ROI is useful to
%         exclude partial volume effects of voxels the surfaces pass
%         through, and will cause the mapping to ignore voxels that don't
%         have a positive value in the mask.  The subdivision number
%         specifies how it approximates the amount of the volume the
%         polyhedron intersects, by splitting each voxel into NxNxN pieces,
%         and checking whether the center of each piece is inside the
%         polyhedron.  If you have very large voxels, consider increasing
%         this if you get zeros in your output.
%
%
%
% Written by E. Gordon 2/15/13

%directory the surfaces are in
surfdir = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/';

%where workbench is installed
workbenchdir = '';%'/data/cn4/evan/workbench/bin_linux64/';

%where the 711-2B to MNI tranform is located on the server
T4_7112B_to_MNI = '/data/nil-bluearc/raichle/lin64-tools/711-2B_to_MNI152lin_T1_t4';

%paths to the left and right surfaces as well as the white and pial surfaces
Lsurface = [surfdir 'Conte69.L.midthickness.32k_fs_LR.surf.gii'];
Lwhitesurface = [surfdir 'Conte69.L.white.32k_fs_LR.surf.gii'];
Lpialsurface = [surfdir 'Conte69.L.pial.32k_fs_LR.surf.gii'];
Rsurface = [surfdir 'Conte69.R.midthickness.32k_fs_LR.surf.gii'];
Rwhitesurface = [surfdir 'Conte69.R.white.32k_fs_LR.surf.gii'];
Rpialsurface = [surfdir 'Conte69.R.pial.32k_fs_LR.surf.gii'];




%Set defaults

if ~exist('mappingtype') || isempty(mappingtype)
    mappingtype = 'ribbon-constrained';
end

if ~exist('hem') || isempty(hem)
    hem = 'both';
end

if ~exist('space')
    space = 'MNI';
end

mappingtype_orig = mappingtype;

thisdir = pwd;

volumes = dir(volumepattern);

for volumenum = 1:length(volumes)
    
    slashlocs = strfind(volumepattern,'/');
    if isempty(slashlocs)
        volume = [thisdir '/' volumes(volumenum).name];
    else
        volume = [volumepattern(1:slashlocs(end)) volumes(volumenum).name];
    end
    
    deletenifti = 0;
    deleteMNInifti = 0;
    
    %check if volume is a 4dfp and convert to nifti if needed
    if strcmp(volume(end-7:end),'4dfp.img')
        system(['nifti_4dfp -n ' volume ' ' volume(1:end-8) 'nii'])
        volume = [volume(1:end-8) 'nii'];
        deletenifti = 1;
    elseif strcmp(volume(end-6:end),'.nii.gz')
        gunzip(volume)
        volume = [volume(1:end-3)];
        deletenifti = 1;
    end
    
    %Transform to MNI if desired
    if strcmp(space,'711-2B')
        system(['nifti_4dfp -4 ' volume ' Temp.4dfp.img'])
        system(['t4img_4dfp ' T4_7112B_to_MNI ' Temp.4dfp.img Temp_MNI.4dfp.img'])
        system(['nifti_4dfp -n Temp_MNI.4dfp.img ' volume(1:end-4) '_MNI.nii'])
        delete('Temp*.4dfp*')
        deleteMNInifti = 2;
        volume = [volume(1:end-4) '_MNI.nii'];
    end
    
    
    %Map to surface
        
    if strcmp(hem,'both') || strcmp(hem,'L')
        
        if strcmp(mappingtype_orig,'ribbon-constrained')
            mappingtype = ['ribbon-constrained ' Lwhitesurface ' ' Lpialsurface ' -voxel-subdiv 5'];
        end
        
        system([workbenchdir 'wb_command -volume-to-surface-mapping ' volume ' ' Lsurface ' ' volume(1:end-4) '_L.func.gii -' mappingtype]);
        
    end
    
    if strcmp(hem,'both') || strcmp(hem,'R')
        
        if strcmp(mappingtype_orig,'ribbon-constrained')
            mappingtype = ['ribbon-constrained ' Rwhitesurface ' ' Rpialsurface ' -voxel-subdiv 5'];
        end
        
        system([workbenchdir 'wb_command -volume-to-surface-mapping ' volume ' ' Rsurface ' ' volume(1:end-4) '_R.func.gii -' mappingtype]);
        
    end
    
    
    %Clean up temporary files
    
    if deletenifti
        delete(volume)
    end
    if deleteMNInifti
        warning off
        delete(volume)
        if deletenifti
            delete([volume(1:end-8) '.nii'])
        end
    end
    
end


