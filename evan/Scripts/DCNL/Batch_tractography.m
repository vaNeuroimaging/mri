%Batch_tractography
%
%For any number of subjects, conducts tractography between any number of
%pairs of seed images. Note that this is pairwise seed-to-seed
%tractography (fancy stuff like waypoints not supported). Resulting tracts
%will be in subject (diffusion) space rather than standard space.
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
subjects = {'208','222','227','253','256','258','270','292','301','305','309','343','362','374','383','395','415','416','417','420'};
%{'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%;
%{'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};

%,'309'

%'112','126','133','137','208','222','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420',
%

%List of seeds with full pathnames
seeds = {'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/DC_2mm.nii',...
'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Chelsea/NewTLT_caudate_vs_Overall_IALscore/Hip_cluster_TLT.nii',...
'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_PFC.nii'};


%{'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/DC_2mm.nii',...
 %   '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest_old/2mm/dmPFC.nii',...
  %  '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest_old/2mm/laIns.nii',...
   % '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest_old/2mm/ldlPFC.nii',...
    %'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest_old/2mm/lPar.nii',...
    %'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest_old/2mm/raIns.nii',...
    %'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest_old/2mm/rdlPFC.nii',...
    %'/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest_old/2mm/avoid_brainstem.nii'};

%Pairwise tractography to run.  Each bracketed number pair is tractography
%between one pair of seeds; each number within the brackets is a seed (the
%number reflects the order within the "seeds" variable)
seedconnections = {[1 2]};


%Paths threshold of fdt image. All voxels with path counts less than
%(pathsthreshold * the max path count) will be masked out
pathsthreshold = .15;

%Location where mean FA values within tracts will be written out
outputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Chelsea/MeanFAOutput.txt';


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
    eval(['!flirt -in ' basedir 'bothruns/nodif_brain.nii.gz -ref /data/apps/fsl/4.1.9/data/standard/MNI152_T1_2mm_brain -omat ' basedir 'bothruns.bedpostX/xfms/diff2standard.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
    
    %invert that to get the standard space to diffusion space transformation
    eval(['!convert_xfm -omat ' basedir 'bothruns.bedpostX/xfms/standard2diff.mat -inverse ' basedir 'bothruns.bedpostX/xfms/diff2standard.mat']);
    
    %determine diffusion space to the standard space transformation
    eval(['!flirt -in ' basedir 'bothruns/nodif_brain.nii.gz -ref ' basedir '../SPM8/FirstRest/fsw0000.nii -omat ' basedir 'bothruns.bedpostX/xfms/diff2func.mat -bins 256 -cost normmi -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -dof 12  -interp trilinear']);
    
    %invert that to get the standard space to diffusion space transformation
    eval(['!convert_xfm -omat ' basedir 'bothruns.bedpostX/xfms/func2diff.mat -inverse ' basedir 'bothruns.bedpostX/xfms/diff2func.mat']);
    
    for seednum = 1:length(seeds)
    
        %Get the names of the seed images
        thisseedname = seeds{seednum};
        slashpositions = find(thisseedname=='/');
        shortseedname{seednum} = thisseedname(slashpositions(end)+1:end-4);
        transformedseedname{seednum} = [basedir 'diffspace_' shortseedname{seednum} '.nii.gz'];
        
        %Transform the seed images into diffusion space
        if seednum == 2
            
            if exist([basedir '/func2diff.mat'], 'file')
                
                eval(['!flirt -in ' thisseedname ' -ref ' basedir 'bothruns/nodif_brain.nii.gz  -applyxfm -init ' basedir '/func2diff.mat -out ' transformedseedname{seednum} ' -interp nearestneighbour']);
                
            else
                eval(['!flirt -in ' thisseedname ' -ref ' basedir 'bothruns/nodif_brain.nii.gz  -applyxfm -init ' basedir 'bothruns.bedpostX/xfms/func2diff.mat -out ' transformedseedname{seednum} ' -interp nearestneighbour']);
            end
            
        else
            
            eval(['!flirt -in ' thisseedname ' -ref ' basedir 'bothruns/nodif_brain.nii.gz  -applyxfm -init ' basedir 'bothruns.bedpostX/xfms/standard2diff.mat -out ' transformedseedname{seednum} ' -interp nearestneighbour']);
        end
        
%         %Threshold diffusion-space masks at .5 (to deal with interpolation) and binarize
%         eval(['!fslmaths ' transformedseedname{seednum} ' -thr .5 -bin ' transformedseedname{seednum}]);
        
    end
    
    for seedconnection = 1:length(seedconnections)
        %disp([shortseedname{seedconnections{seedconnection}(1)} ' to ' shortseedname{seedconnections{seedconnection}(2)} ' via ' shortseedname{3}])
        disp([shortseedname{seedconnections{seedconnection}(1)} ' to ' shortseedname{seedconnections{seedconnection}(2)}])
        
        %Make the output folder
        %outfolder = [basedir 'bothruns_' shortseedname{seedconnections{seedconnection}(1)} '_to_' shortseedname{seedconnections{seedconnection}(2)} '_via_' shortseedname{3} '/'];
        outfolder = [basedir 'bothruns_' shortseedname{seedconnections{seedconnection}(1)} '_to_' shortseedname{seedconnections{seedconnection}(2)} '/'];
        try rmdir(outfolder,'s');catch;end
        mkdir(outfolder);
    
        %Run tractography
        temp = evalc(['!/data/apps/fsl/4.1.9/bin/probtrackx --network --mode=seedmask -x ' transformedseedname{seedconnections{seedconnection}(2)} ' --waypoints=' transformedseedname{seedconnections{seedconnection}(1)} ' --stop=' transformedseedname{seedconnections{seedconnection}(1)} ' -l -c 0.1 -S 2000 --steplength=0.5 -P 5000 --forcedir --opd -s ' basedir 'bothruns.bedpostX/merged -m ' basedir 'bothruns.bedpostX/nodif_brain_mask  --dir=' outfolder]);  %' --avoid=' transformedseedname{end}
        %temp = evalc(['!/data/apps/fsl/4.1.9/bin/probtrackx --network --mode=seedmask -x ' seeds{seedconnections{seedconnection}(2)} ' --waypoints=' seeds{seedconnections{seedconnection}(1)} ' --stop=' seeds{seedconnections{seedconnection}(1)} ' --avoid=' seeds{end} ' -l -c 0.1 -S 200 --xfm=' basedir 'bothruns.bedpostX/xfms/standard2diff.mat --steplength=0.5 -P 10000 --forcedir --opd -s ' basedir 'bothruns.bedpostX/merged -m ' basedir 'bothruns.bedpostX/nodif_brain_mask  --dir=' outfolder]);
        
        %Find the maximum number of paths that passed through any one voxel
        pathcountoutput = evalc(['!fslstats ' outfolder 'fdt_paths.nii.gz -R']);
        minmaxpaths = str2num(pathcountoutput);
        maxpaths = minmaxpaths(2);
        
        %Threshold the tract by a percentage (defined by pathsthreshold) of the maximum path number, and binarize
        eval(['!fslmaths ' outfolder 'fdt_paths.nii.gz -thr ' num2str(maxpaths*pathsthreshold) ' -bin ' outfolder 'fdt_paths_thresh_mask']);
        
        %Calculate the volume of the thresholded tract
        voxelsvolumeoutput = evalc(['!fslstats ' outfolder 'fdt_paths_thresh_mask -V']);
        voxelsvolume = str2num(voxelsvolumeoutput);
        volume = num2str(voxelsvolume(2));
        
        %Mask the subject's FA image by the binarized tract
        eval(['!fslmaths ' basedir 'bothruns/dti_FA.nii.gz -mas ' outfolder 'fdt_paths_thresh_mask.nii.gz ' outfolder 'FA_within_tract']);
        
        %Caluclate the mean FA within this tract-masked FA image
        meanFA = evalc(['!fslstats ' outfolder 'FA_within_tract.nii.gz -M']);
        
        %Write out mean FA and volume data for this tract
        %texttowrite = [subjid,'   ', [shortseedname{seedconnections{seedconnection}(1)} '_to_' shortseedname{seedconnections{seedconnection}(2)} '_via_' shortseedname{3}],'   ',num2str(str2num(meanFA)),'   ',num2str(volume)];
        texttowrite = [subjid,'   ', [shortseedname{seedconnections{seedconnection}(1)} '_to_' shortseedname{seedconnections{seedconnection}(2)}],'   ',num2str(str2num(meanFA)),'   ',num2str(volume)];
        dlmwrite([outputfilename],texttowrite,'-append','delimiter','');

    end
    
end