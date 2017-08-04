subjects = {'396'};

seednames = {'PCC'};

seedimages = {'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/PCC_aal_funcspace.nii'};

for sub = 1:length(subjects)
    subjid = subjects{sub};


fiberfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/Fiber.dat'];


data_location = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/SPM8/Rest/'];
datafile4D =  [data_location 'filt_smoothed.nii.gz'];
datafiles3D = [data_location 'fsw*.nii'];
%normdatafiles3D = [data_location 'waD*.hdr'];
outputfile = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/'];

DTIdir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/'];


datafilenames3D = dir(datafiles3D);
m = size(datafilenames3D, 1);
for i=1:m
    P(i, :) = [data_location datafilenames3D(i).name];
end

% normdatafilenames3D = dir(normdatafiles3D);
% m = size(normdatafilenames3D, 1);
% for i=1:m
%     nP(i, :) = [data_location normdatafilenames3D(i).name];
% end


for seednum = 1:length(seednames)
    seedimage = seedimages{seednum};
    seedname = seednames{seednum};
    
    
    
    %VOXELWISE STRUCTURAL CONNECTIVITY
    %----------------------------------------------------------------------
    
    
    
    eval(['!rm ' DTIdir 'bothruns/' seedname '_diffspace*']);
eval(['!flirt -in ' P(1, :) ' -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -omat ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
eval(['!convert_xfm -omat ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -inverse ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat']);
eval(['!flirt -ref ' DTIdir 'bothruns/dti_FA.nii.gz -in ' seedimage ' -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -out ' DTIdir 'bothruns/' seedname '_diffspace']);
eval(['!fslmaths ' DTIdir 'bothruns/' seedname '_diffspace.nii.gz -thr .5 -bin ' DTIdir 'bothruns/' seedname '_diffspace.nii.gz']);
gunzip([DTIdir 'bothruns/' seedname '_diffspace.nii.gz']);


loadedseed = load_nii([DTIdir 'bothruns/' seedname '_diffspace.nii']);
seeddata = loadedseed.img;

fulldataimage = cell(size(seeddata));


fid = fopen(fiberfilename);

%Header
nothing = fread(fid,8,'schar=>char')';
nFibers = fread(fid,1,'int32');
nothing = fread(fid,1,'int32');
nothing = fread(fid,1,'float32');
imgx = fread(fid,1,'int32');
imgy = fread(fid,1,'int32');
imgz = fread(fid,1,'int32');
voxelsizex = fread(fid,1,'float32');
voxelsizey = fread(fid,1,'float32');
voxelsizez = fread(fid,1,'float32');

fseek(fid,128,-1);
% test = zeros(10000000,3);
% testcounter = 1;

for fiber = 1:nFibers
    fibers(fiber).length = fread(fid,1,'int32');
    nothing = fread(fid,1,'uchar=>char');
    nothing = fread(fid,1,'uint8');
    nothing = fread(fid,1,'uint8');
    nothing = fread(fid,1,'uint8');
    
    fibers(fiber).start = fread(fid,1,'int32');
    fibers(fiber).end =fread(fid,1,'int32');
    fiberpointcounter = 1;
    for fiberpoint = 1:fibers(fiber).length
        x = floor(fread(fid,1,'float32')+1);
        y = floor(fread(fid,1,'float32')+1);
        z = floor(fread(fid,1,'float32')+1);
        if x>0 && y>0 && z>0 && x<=imgx && y<=imgy && z<=imgz
            fibers(fiber).points(fiberpointcounter,:) = [x y z];
            fulldataimage{x,y,z} = [fulldataimage{x,y,z} fiber];
            fiberpointcounter = fiberpointcounter+1;
        end
        %test(testcounter,:) = [x y z];
        %testcounter = testcounter+1;
        
    end
end

fclose(fid);

seedvoxels = [];
for z = 1:size(seeddata,3)
    [x, y] = find(squeeze(seeddata(:,:,z)));
    seedvoxels = [seedvoxels; [x y repmat(z,length(x),1)]];
end

fibersthroughseed = [];
for seedvoxel = 1:size(seedvoxels,1)
    fibersthroughseed = [fibersthroughseed fulldataimage{seedvoxels(seedvoxel,1),seedvoxels(seedvoxel,2),seedvoxels(seedvoxel,3)}];
end
fibersthroughseed = unique(fibersthroughseed);




%FIBER DENSITY

fiberdensityimage = zeros(size(seeddata));
for fibercount = 1:length(fibersthroughseed)
    fibernum = fibersthroughseed(fibercount);
    for pointcount = 1:fibers(fibernum).length
        coords = fibers(fibernum).points(pointcount,:);
        fiberdensityimage(coords(1),coords(2),coords(3)) = fiberdensityimage(coords(1),coords(2),coords(3))+1;
    end
end

%Make the output nifti
fiberdensityoutput = make_nii(fiberdensityimage, loadedseed.hdr.dime.pixdim(2:4), loadedseed.hdr.hist.originator(1:3));

%Save the output nifti
save_nii(fiberdensityoutput,[outputfile 'Fiberdensity_to_' seedname '.nii']);





%FA

gunzip([DTIdir 'bothruns/dti_FA.nii.gz']);
FAimagefile = load_nii([DTIdir 'bothruns/dti_FA.nii']);
FAimage = FAimagefile.img;
clear FAimagefile

FAoutputimage = zeros(size(seeddata));

for x = 1:size(seeddata,1)
    for y = 1:size(seeddata,2)
        for z= 1:size(seeddata,3)
            
            if fiberdensityimage(x,y,z)>0 && seeddata(x,y,z)==0;
                
                allpointsbetween = [];
                
                for fibercount = 1:length(fulldataimage{x,y,z})
                    thisfiber = fulldataimage{x,y,z}(fibercount);
                    
                    seedvoxellocations = [];
                    for seedvoxel = 1:size(seedvoxels,1)
                        seedvoxellocations = [seedvoxellocations ...
                            intersect(find(fibers(thisfiber).points(:,1)==seedvoxels(seedvoxel,1)),find(fibers(thisfiber).points(:,2)==seedvoxels(seedvoxel,2)),find(fibers(thisfiber).points(:,3)==seedvoxels(seedvoxel,3)))];
                    end
                    
                    targetvoxellocation = intersect(find(fibers(thisfiber).points(:,1)==x),find(fibers(thisfiber).points(:,2)==y),find(fibers(thisfiber).points(:,3)==z));
                    
                    [mindistance closestvoxel] = min(abs(seedvoxellocations-targetvoxellocation));
                    if targetvoxellocation-closestvoxel<0
                        fiberpointsbetween = fibers(thisfiber).points(targetvoxellocation:closestvoxel,:);
                    else
                        fiberpointsbetween = fibers(thisfiber).points(closestvoxel:targetvoxellocation,:);
                    end
                    
                    allpointsbetween = [allpointsbetween; fiberpointsbetween];
                end
                
                allpointsbetween = unique(allpointsbetween,'rows');
                
                FAvals = [];
                for point = 1:size(allpointsbetween)
                    FAvals = [FAvals FAimage(allpointsbetween(point))];
                end
                
                FAoutputimage(x,y,z) = mean(FAvals);
            end
        end
    end
end
                        
%Make the output nifti
FAoutput = make_nii(FAoutputimage, loadedseed.hdr.dime.pixdim(2:4), loadedseed.hdr.hist.originator(1:3));

%Save the output nifti
save_nii(FAoutput,[outputfile 'FA_to_' seedname '.nii']);     

end

end
                        
