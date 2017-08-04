%Batch_tractography
%
%For any number of subjects, conducts tractography from any number of
% seed images. Resulting tracts will be in subject (diffusion) space rather
% than standard space.
%
%The script will then threshold the resulting tract, calculate the volume
%of the tract and the mean FA value within that tract, and save those
%values to a text file which can be easily imported into Excel.
%
%Subjects/seeds are specified at the top of the script. All subjects' DTI
%data must be fully preprocessed, up through bedposting. All seeds must be
%in the same space as the MNI_152_2mm_brain.nii.gz standard image provided in FSL.
%
%You must have FSL loaded on your computer for this to run.
%
%Created by E. Gordon 04/11



%USER INPUT (and more down below)
%--------------------------------------------------------------------------

%List of subjects
subjects = {'126','112','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
    %'113','101','102','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};

%List of seeds with full pathnames
seeds = {'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/DC_2mm.nii'};


%Paths threshold of fdt image. All voxels with path counts less than
%(pathsthreshold * the max path count) will be masked out
pathsthreshold = .1;

%Location where mean FA values within tracts will be written out
outputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Tractography/MeanFAOutput_oneseed.txt';


%END USER INPUT 
%--------------------------------------------------------------------------


delete(outputfilename);
fid = fopen(outputfilename,'at');
fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Subject','Tract','meanFA','Volume');
fclose(fid);
dlmwrite([ outputfilename],' ','-append');


for subject = 1:length(subjects)
    subjid = subjects{subject};
    disp(subjid)
    
    
    %USER INPUT
    %--------------------------------------------------------------------------
    
    %DTI directory for this subject
    basedir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/'];
    
    %END USER INPUT
    %--------------------------------------------------------------------------
    
    %determine diffusion space to the standard space transformation
    eval(['!flirt -in ' basedir 'bothruns/nodif_brain.nii.gz -ref /data/apps/fsl/4.1.7/data/standard/MNI152_T1_2mm_brain -omat ' basedir 'bothruns.bedpostX/xfms/diff2standard.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
    
    %invert that to get the standard space to diffusion space transformation
    eval(['!convert_xfm -omat ' basedir 'bothruns.bedpostX/xfms/standard2diff.mat -inverse ' basedir 'bothruns.bedpostX/xfms/diff2standard.mat']);
    
    for seednum = 1:length(seeds)
    
        %Get the names of the seed images
        thisseedname = seeds{seednum};
        slashpositions = find(thisseedname=='/');
        shortseedname{seednum} = thisseedname(slashpositions(end)+1:end-4);
        transformedseedname{seednum} = [basedir 'diffspace_' shortseedname{seednum} '.nii.gz'];
        
        %Transform the seed images into diffusion space
        eval(['!flirt -in ' thisseedname ' -ref ' basedir 'bothruns/nodif_brain.nii.gz  -applyxfm -init ' basedir 'bothruns.bedpostX/xfms/standard2diff.mat -out ' transformedseedname{seednum}]);
        
        %Threshold diffusion-space masks at .5 (to deal with interpolation) and binarize
        eval(['!fslmaths ' transformedseedname{seednum} ' -thr .5 -bin ' transformedseedname{seednum}]);
        
    
        disp([shortseedname{seednum}])
        
        %Make the output folder
        outfolder = [basedir 'bothruns_' shortseedname{seednum} '/'];
        try rmdir(outfolder,'s');catch;end
        mkdir(outfolder);
    
%         %Make the masks.txt file (a list of masks) used by the tractography program
%         delete([outfolder 'masks.txt']);
%         fid = fopen([outfolder 'masks.txt'],'at');
%         fprintf(fid,'%s\n\r\%s',transformedseedname{1},transformedseedname{2});
%         fclose(fid);
    
        %Run tractography
        temp = evalc(['!/data/apps/fsl/4.1.7/bin/probtrackx --network --mode=seedmask -x ' transformedseedname{seednum} ' -l -c 0.2 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd -s ' basedir 'bothruns.bedpostX/merged -m ' basedir 'bothruns.bedpostX/nodif_brain_mask  --dir=' outfolder]);
        
        %Find the maximum number of paths that passed through any one voxel
        pathcountoutput = evalc(['!fslstats ' outfolder 'fdt_paths.nii.gz -R']);
        minmaxpaths = str2num(pathcountoutput);
        maxpaths = minmaxpaths(2);
        
        %Threshold the tract by a percentage (defined by pathsthreshold) of the maximum path number, and binarize
        eval(['!fslmaths ' outfolder 'fdt_paths.nii.gz -thr ' num2str(maxpaths*pathsthreshold) ' -bin ' outfolder 'fdt_paths_thresh_mask']);
        
        %Calculate the volume of the thresholded tract
        voxelsvolumeoutput = evalc(['!fslstats ' outfolder 'fdt_paths_thresh_mask -V']);
        voxelsvolume = str2num(voxelsvolumeoutput);
        volume = voxelsvolume(2);
        
        %Mask the subject's FA image by the binarized tract
        eval(['!fslmaths ' basedir 'bothruns/dti_FA.nii.gz -mas ' outfolder 'fdt_paths_thresh_mask.nii.gz ' outfolder 'FA_within_tract']);
        
        %Caluclate the mean FA within this tract-masked FA image
        meanFA = evalc(['!fslstats ' outfolder 'FA_within_tract.nii.gz -M']);
        
        %Write out mean FA and volume data for this tract
        texttowrite = [subjid,'   ', [shortseedname{seednum}],'   ',num2str(meanFA),'   ',num2str(volume)];
        dlmwrite([outputfilename],texttowrite,'-append','delimiter','');

    end
    
end