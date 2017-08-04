
subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','374','383','395','396'};
    %



%'166',

allsubjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
    %


needsaflip = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];

spherevoxelradius = 3;

seednames = {'PCC'};

seedimages = {'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/PCC_aal_funcspace.nii'};

brainmaskname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.hdr';
brainmask = load_nii(brainmaskname);

sublist = str2num(cell2mat(allsubjects'));

warning off

for sub = 1:length(subjects)
    subjid = subjects{sub};


fiberfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/Fiber.dat'];
%FAfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/DTIStudioFA.hdr'];
FAfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/bothruns/dti_FA.nii'];
try gunzip([FAfilename '.gz']); catch; end


data_location = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/SPM8/Rest/'];
datafile4D =  [data_location 'filt_smoothed.nii.gz'];
datafiles3D = [data_location 'fsw*.nii'];
%normdatafiles3D = [data_location 'waD*.hdr'];
outputfile = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/'];

DTIdir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/'];


%     motionparamfile = dir([data_location 'rp*.txt']);
%     motionparams = textread([data_location motionparamfile.name]);


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
    
    FCoutputimgname = [seedname '_FC.nii'];
    
    %VOXELWISE FUNCTIONAL CONNECTIVITY
    %----------------------------------------------------------------------
% 
% 
%     disp(['Subject ' subjid ': Functional Connectivity'])
% 
%     seed = maroi('load_cell',[seedimage(1:end-4) '_roi.mat']);
%     [Y a b c] = getdata(seed{1}, P,'l');
%     seed_timecourse = mean(Y,2);
%     
%     CSF_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjid '_CSF_roi.mat']);
%     [Y a b c] = getdata(CSF_rois{1}, P,'l');
%     CSF_timecourse = mean(Y,2);
%     
%     WM_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjid '_WM_roi.mat']);
%     [Y a b c] = getdata(WM_rois{1}, P,'l');
%     WM_timecourse = mean(Y,2);
%     
%     try gunzip(datafile4D); catch; end;
%     
%     data = load_nii(datafile4D(1:end-3));
%     
%     reshapeddata = double(reshape(data.img,size(data.img,1)*size(data.img,2)*size(data.img,3),size(data.img,4))');
%     
%     %R = jacket_partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams]);
%     R = partialcorr(seed_timecourse, reshapeddata, [CSF_timecourse WM_timecourse motionparams]);
%     Fishervals = .5*(log(1+R)-log(1-R));
%     Fisherimg(:,:,:) = reshape(Fishervals,size(data.img,1),size(data.img,2),size(data.img,3));
%     clear R Fishervals
%     
%     Fisherimg(find(brainmask.img==0)) = 0;
%     
%     try mkdir(outputfile); catch; end
%     
%     try delete([outputfile FCoutputimgname]); catch; end
%     
%     outputimage = make_nii(Fisherimg,brainmask.hdr.dime.pixdim(2:4), brainmask.hdr.hist.originator(1:3));
%     save_nii(outputimage,[outputfile FCoutputimgname]);
%     
%     %----------------------------------------------------------------------
    
    
    
    
    clear Fisherimg data reshapeddata motionparams seed_timecourse WM_timecourse CSF_timecourse Fisheroutput string m
    
    
    
    
    %VOXELWISE STRUCTURAL CONNECTIVITY
    %----------------------------------------------------------------------
    
    disp(['Subject ' subjid ': Structural Connectivity -- reading Fibers file and making Fiber Density map'])
    
evalc(['!rm ' DTIdir seedname '_diffspace*']);
evalc(['!rm ' DTIdir seedname '_flippeddiffspace*']);
evalc(['!rm ' DTIdir 'brainmask_diffspace*']);
evalc(['!rm ' DTIdir 'brainmask_flippeddiffspace*']);


% evalc(['!flirt -in ' DTIdir 'bothruns/dti_FA.nii.gz -ref ' DTIdir 'DTIStudioFA.hdr -omat ' DTIdir 'FSL2DTIS.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear'])
% evalc(['!flirt -in ' DTIdir 'bothruns/nodif_brain.nii.gz -ref ' DTIdir 'DTIStudioFA.hdr -applyxfm -init ' DTIdir 'FSL2DTIS.mat -out ' DTIdir 'nodif_brain_DTISspace'])
% evalc(['!flirt -in ' P(1, :) ' -ref ' DTIdir 'nodif_brain_DTISspace.nii.gz -omat ' DTIdir 'func2DTIS.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear'])
% evalc(['!convert_xfm -omat ' DTIdir 'DTIS2func.mat -inverse ' DTIdir 'func2DTIS.mat'])
% evalc(['!flirt -ref ' DTIdir 'nodif_brain_DTISspace -in ' seedimage ' -applyxfm -init ' DTIdir 'func2DTIS.mat -out ' DTIdir seedname '_diffspace'])

if needsaflip(find(sublist==str2num(subjid)));

eval(['!flirt -in ' P(1, :) ' -ref /fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/' subjid '_nodif_brain.nii.gz -omat ' DTIdir 'bothruns.bedpostX/xfms/func2flippeddiff.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
%eval(['!flirt -in ' P(1, :) ' -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -omat ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);

%eval(['!convert_xfm -omat ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -inverse ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat']);
eval(['!convert_xfm -omat ' DTIdir 'bothruns.bedpostX/xfms/flippeddiff2func.mat -inverse ' DTIdir 'bothruns.bedpostX/xfms/func2flippeddiff.mat']);

eval(['!flirt -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -in ' seedimage ' -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/func2flippeddiff.mat -out ' DTIdir seedname '_flippeddiffspace']);
%eval(['!flirt -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -in ' seedimage ' -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -out ' DTIdir seedname '_diffspace']);
eval(['!flirt -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -in ' brainmaskname ' -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/func2flippeddiff.mat -out ' DTIdir 'brainmask_flippeddiffspace']);


% eval(['!fslmaths ' DTIdir seedname '_diffspace.nii.gz -thr .5 -bin ' DTIdir seedname '_diffspace.nii.gz']);
% gunzip([DTIdir seedname '_diffspace.nii.gz']);
eval(['!fslmaths ' DTIdir seedname '_flippeddiffspace.nii.gz -thr .5 -bin ' DTIdir seedname '_flippeddiffspace.nii.gz']);
gunzip([DTIdir seedname '_flippeddiffspace.nii.gz']);

eval(['!fslmaths ' DTIdir 'brainmask_flippeddiffspace.nii.gz -thr .5 -bin ' DTIdir 'brainmask_flippeddiffspace.nii.gz']);
gunzip([DTIdir 'brainmask_flippeddiffspace.nii.gz']);

loadedseed = load_nii([DTIdir seedname '_flippeddiffspace.nii']);
seeddata = loadedseed.img;

brainmasksub = load_nii([DTIdir 'brainmask_flippeddiffspace.nii']);
brainmaskdatasub = brainmasksub.img;

% loadedaltseed = load_nii([DTIdir seedname '_diffspace.nii']);
% altseeddata = loadedaltseed.img;

else
   
    evalc(['!flirt -in ' P(1, :) ' -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -omat ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
 evalc(['!convert_xfm -omat ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -inverse ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat']);
 evalc(['!flirt -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -in ' seedimage ' -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -out ' DTIdir seedname '_diffspace']);
 eval(['!fslmaths ' DTIdir seedname '_diffspace.nii.gz -thr .5 -bin ' DTIdir seedname '_diffspace.nii.gz']);
gunzip([DTIdir seedname '_diffspace.nii.gz']);
    loadedseed = load_nii([DTIdir seedname '_diffspace.nii']);
seeddata = loadedseed.img;

seeddata = flipdim(flipdim(seeddata,1),2);



evalc(['!flirt -ref ' DTIdir 'bothruns/nodif_brain.nii.gz -in ' brainmaskname ' -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/func2diff.mat -out ' DTIdir 'brainmask_diffspace']);
 eval(['!fslmaths ' DTIdir 'brainmask_diffspace.nii.gz -thr .5 -bin ' DTIdir 'brainmask_diffspace.nii.gz']);
gunzip([DTIdir 'brainmask_diffspace.nii.gz']);
    brainmasksub = load_nii([DTIdir 'brainmask_diffspace.nii']);
brainmaskdatasub = brainmasksub.img;

brainmaskdatasub = flipdim(flipdim(brainmaskdatasub,1),2);

    
end



voxelindex = find(seeddata>-1);
[xindex, yzindex] = find(seeddata>-1);
yindex = rem(yzindex,size(seeddata,2));
yindex(find(yindex==0)) = size(seeddata,2);
zindex = ceil(yzindex/size(seeddata,2));
indexmatrix = uint8([xindex yindex zindex]);



fulldataimage = cell(size(seeddata));

FAimagefile = load_nii(FAfilename);
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
        x = round(fread(fid,1,'float32'))+1;
        y = round(fread(fid,1,'float32'))+1;
        z = round(fread(fid,1,'float32'))+1;
        if x>0 && y>0 && z>0 && x<=imgx && y<=imgy && z<=imgz
            %fibers(fiber).points(fiberpointcounter,:) = [x y z];
            fibers(fiber).points(fiberpointcounter) = (z-1)*size(seeddata,1)*size(seeddata,2) + (y-1)*size(seeddata,1) + x;
            fulldataimage{x,y,z} = [fulldataimage{x,y,z} fiber];
            fiberpointcounter = fiberpointcounter+1;
        end
        %test(testcounter,:) = [x y z];
        %testcounter = testcounter+1;
        
    end
end

fclose(fid);

% seedvoxels = [];
% for z = 1:size(seeddata,3)
%     [x, y] = find(squeeze(seeddata(:,:,z)));
%     seedvoxels = [seedvoxels; [x y repmat(z,length(x),1)]];
% end
seedvoxels = find(seeddata);

fibersthroughseed = [];
for seedvoxel = 1:length(seedvoxels)
    fibersthroughseed = [fibersthroughseed fulldataimage{indexmatrix(seedvoxels(seedvoxel),1),indexmatrix(seedvoxels(seedvoxel),2),indexmatrix(seedvoxels(seedvoxel),3)}];
end
fibersthroughseed = unique(fibersthroughseed);




%FIBER DENSITY
fiberpresencematrix = zeros(size(indexmatrix,1),500);
fiberdensityimage = zeros(size(seeddata));
for fibercount = 1:length(fibersthroughseed)
    fibernum = fibersthroughseed(fibercount);
    for pointnum = 1:length(fibers(fibernum).points)
        fiberdensityimage(fibers(fibernum).points(pointnum)) = fiberdensityimage(fibers(fibernum).points(pointnum))+1;
        previousfibershere = length(find(fiberpresencematrix(fibers(fibernum).points(pointnum),:)));
        fiberpresencematrix(fibers(fibernum).points(pointnum),previousfibershere+1) = fibernum;
    end
end

if ~needsaflip(find(sublist==str2num(subjid)));
    fiberdensityoutputimage = flipdim(flipdim(fiberdensityimage,1),2);
else
    fiberdensityoutputimage = fiberdensityimage;
end

%Make the output nifti
fiberdensityoutput = make_nii(fiberdensityoutputimage, loadedseed.hdr.dime.pixdim(2:4), loadedseed.hdr.hist.originator(1:3));

evalc(['!rm ' outputfile 'Fiberdensity_to_' seedname '*']);

%Save the output nifti
save_nii(fiberdensityoutput,[outputfile 'Fiberdensity_to_' seedname '_diffspace.nii']);

eval(['!flirt -ref ' P(1, :) ' -in ' outputfile 'Fiberdensity_to_' seedname '_diffspace.nii -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -out ' outputfile 'Fiberdensity_to_' seedname '.nii.gz']);
gunzip([outputfile 'Fiberdensity_to_' seedname '.nii.gz']);





% %FA
disp(['Subject ' subjid ': Structural Connectivity -- making voxelwise FA map'])
%gunzip([DTIdir 'bothruns/dti_FA.nii.gz']);

FAimage = FAimagefile.img;



%if ~needsaflip(find(sublist==str2num(subjid)));
    FAimage = flipdim(flipdim(FAimage,1),2);
%end

clear FAimagefile

FAoutputimage = zeros(size(seeddata));
%pointsbetweenimage = cell(size(seeddata));
pointsbetweenmatrix = zeros(size(indexmatrix,1),500);

for voxelnum = 1:size(indexmatrix,1)
    x = indexmatrix(voxelnum,1);
    y = indexmatrix(voxelnum,2);
    z = indexmatrix(voxelnum,3);

    
    
    
    if fiberdensityimage(x,y,z)>0 && seeddata(x,y,z)==0;
        %disp([x y z])
        allpointsbetween = [];
        
        for fibercount = 1:length(fulldataimage{x,y,z})
            thisfiber = fulldataimage{x,y,z}(fibercount);
            if ~isempty(find(fibersthroughseed==thisfiber,1))
                
                
%                 seedvoxellocations = [];
%                 for seedvoxel = 1:size(seedvoxels,1)
%                     
%                     seedvoxellocations = [seedvoxellocations; ...
%                     intersect(intersect(find(fibers(thisfiber).points(:,1)==seedvoxels(seedvoxel,1)),find(fibers(thisfiber).points(:,2)==seedvoxels(seedvoxel,2))),find(fibers(thisfiber).points(:,3)==seedvoxels(seedvoxel,3)))];
%                 end
                
                seedvoxellocations = find(ismember(fibers(thisfiber).points,seedvoxels));
                
                
%                targetvoxellocation = intersect(intersect(find(fibers(thisfiber).points(:,1)==x),find(fibers(thisfiber).points(:,2)==y)),find(fibers(thisfiber).points(:,3)==z));
                targetvoxellocation = find(ismember(fibers(thisfiber).points,voxelnum));
                
                [mindistance closestvoxelindex] = min(abs(seedvoxellocations-targetvoxellocation(1)));
                closestvoxel = seedvoxellocations(closestvoxelindex);
                if targetvoxellocation-closestvoxel<0
                    fiberpointsbetween = fibers(thisfiber).points(targetvoxellocation:closestvoxel);
                else
                    fiberpointsbetween = fibers(thisfiber).points(closestvoxel:targetvoxellocation);
                end
                
                
                fiberpointindices = find(ismember((1:size(indexmatrix,1))',fiberpointsbetween));
                
                allpointsbetween = [allpointsbetween; fiberpointindices];
            end
        end
        
        allpointsbetween = unique(allpointsbetween);
        
        pointsbetweenmatrix(voxelnum,1:length(allpointsbetween)) = allpointsbetween';
        
        FAvals = FAimage(allpointsbetween);
        
%         FAvals = [];
%         for point = 1:length(allpointsbetween)
%             FAvals = [FAvals FAimage(allpointsbetween(point))];
%         end
        
        FAoutputimage(x,y,z) = mean(FAvals);
    end
    
end

if ~needsaflip(find(sublist==str2num(subjid)));
    FAoutputimage = flipdim(flipdim(FAoutputimage,1),2);
end
                        
%Make the output nifti
FAoutput = make_nii(FAoutputimage, loadedseed.hdr.dime.pixdim(2:4), loadedseed.hdr.hist.originator(1:3));

evalc(['!rm ' outputfile 'FA_to_' seedname '*']);

%Save the output nifti
save_nii(FAoutput,[outputfile 'FA_to_' seedname '_diffspace.nii']);

eval(['!flirt -ref ' P(1, :) ' -in ' outputfile 'FA_to_' seedname '_diffspace.nii -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -out ' outputfile 'FA_to_' seedname '.nii.gz']);
eval(['!fslmaths ' outputfile 'FA_to_' seedname '.nii.gz -thr .1 ' outputfile 'FA_to_' seedname '.nii.gz']);
gunzip([outputfile 'FA_to_' seedname '.nii.gz']);




disp(['Subject ' subjid ': Structural Connectivity -- making spherewise FA and Fiber Density maps'])



gindexmatrix = gint8(indexmatrix);
FAoutputimagesphere = zeros(size(seeddata));
Fiberdensityoutputimagesphere = zeros(size(seeddata));
%gpointsbetweenmatrix = guint32(pointsbetweenmatrix);
spherevoxelradius = gsingle(spherevoxelradius);

for x = 1:size(seeddata,1)
    for y = 1:size(seeddata,2)
        for z= 1:size(seeddata,3)
        
            if brainmaskdatasub(x,y,z)
        
        distances = gsingle([gindexmatrix(:,1)-x gindexmatrix(:,2)-y gindexmatrix(:,3)-z]);
        euclidean = sqrt(sum((distances.^2),2));
        pointswithinsphere = double(find(euclidean<=spherevoxelradius));
        
        fibersthroughsphere = fiberpresencematrix(pointswithinsphere,:);
        if ~isempty(find(fibersthroughsphere,1))
            Fiberdensityoutputimagesphere(x,y,z) = length(unique(fibersthroughsphere(find(fibersthroughsphere))));
        end
        
        
        pointsbetween2D = pointsbetweenmatrix(pointswithinsphere,:);
        
        if ~isempty(find(pointsbetween2D,1))
        
        allpointsbetween = reshape(pointsbetween2D,(size(pointsbetween2D,1)*size(pointsbetween2D,2)),1);
        allpointsbetween = unique(allpointsbetween);
        allpointsbetween = allpointsbetween(find(allpointsbetween));
        
        
        FAvals = FAimage(allpointsbetween);
        
        FAoutputimagesphere(x,y,z) = mean(FAvals);
        end
            end
        end
    end
    
end
            
            
            

if ~needsaflip(find(sublist==str2num(subjid)));
    FAoutputimagesphere = flipdim(flipdim(FAoutputimagesphere,1),2);
    Fiberdensityoutputimagesphere = flipdim(flipdim(Fiberdensityoutputimagesphere,1),2);
end
                        
%Make the output nifti
FAoutput = make_nii(FAoutputimagesphere, loadedseed.hdr.dime.pixdim(2:4), loadedseed.hdr.hist.originator(1:3));
Fiberdensityoutput = make_nii(Fiberdensityoutputimagesphere, loadedseed.hdr.dime.pixdim(2:4), loadedseed.hdr.hist.originator(1:3));



%Save the output niftis
save_nii(FAoutput,[outputfile 'FA_to_' seedname '_diffspace_sphere.nii']);

eval(['!flirt -ref ' P(1, :) ' -in ' outputfile 'FA_to_' seedname '_diffspace_sphere.nii -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -out ' outputfile 'FA_to_' seedname '_sphere.nii.gz']);
eval(['!fslmaths ' outputfile 'FA_to_' seedname '_sphere.nii.gz -thr .1 ' outputfile 'FA_to_' seedname '_sphere.nii.gz']);
gunzip([outputfile 'FA_to_' seedname '_sphere.nii.gz']);


save_nii(Fiberdensityoutput,[outputfile 'Fiberdensity_to_' seedname '_diffspace_sphere.nii']);

eval(['!flirt -ref ' P(1, :) ' -in ' outputfile 'Fiberdensity_to_' seedname '_diffspace_sphere.nii -applyxfm -init ' DTIdir 'bothruns.bedpostX/xfms/diff2func.mat -out ' outputfile 'Fiberdensity_to_' seedname '_sphere.nii.gz']);
eval(['!fslmaths ' outputfile 'Fiberdensity_to_' seedname '_sphere.nii.gz -thr 1 ' outputfile 'Fiberdensity_to_' seedname '_sphere.nii.gz']);
gunzip([outputfile 'Fiberdensity_to_' seedname '_sphere.nii.gz']);

%----------------------------------------------------------------------


clear fibers FAimage fibersthroughseed  seedvoxels fulldataimage fiberdensityimage seedvoxellocations allpointsbetween FAvals FAoutputimage FAimage

end

end
                        
