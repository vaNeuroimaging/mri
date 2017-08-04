warning off

addpath /data/cn/data1/scripts/fcimage_analysis/fcimage_analysis_v2

subjects = {'vc34096','vc33457' 'vc34125' 'vc34126' 'vc34128' 'vc34140' 'vc34141' 'vc34198' 'vc34199' 'vc34200' 'vc34201' 'vc34220' 'vc33378' 'vc35175' 'vc34252' 'vc33775_2' 'vc35469' 'vc34306' 'vc34307' 'vc34308' 'vc34330' 'vc34331' 'vc34401' 'vc34402' 'vc34403' 'vc34404' 'vc33769' 'vc34408'};
% 

seeds = {'ArticPoint_dmPFC','ArticPoint_rPremot','ArticPoint_raI','Twonetworkpoint_mPFC','Twonetworkpoint_rPar','Twonetworkpoint_rSM','Onenetworkpoint_Prec','Onenetworkpoint_Vis','Onenetworkpoint_dmPFC'};

seedfolder = '/data/cn4/evan/ROIs/';

windowsize = 30;

voxelspace = '333';

brainmaskname = ['/home/usr/fidl/lib/glm_atlas_mask_' voxelspace '.4dfp.img'];

brainmask = read_4dfpimg(brainmaskname);
brainmask = logical(brainmask);

networkimagename = '/data/cn4/evan/ROIs/Power_et_al_2011_Neuron_711.4dfp.img';

networks = read_4dfpimg_HCP(networkimagename);

networks = networks(brainmask);

numnetworks = max(max(networks));


for seednum = 1:length(seeds)

        try seedmask{seednum} = read_4dfpimg_HCP([seedfolder seeds{seednum} '.4dfp.img']);
        catch
            evalc(['!nifti_4dfp -4 ' seedfolder seeds{seednum} '.nii ' seedfolder seeds{seednum} '.4dfp.img']);
            seedmask{seednum} = read_4dfpimg_HCP([seedfolder seeds{seednum} '.4dfp.img']);
        end
        
        seedmask{seednum} = logical(seedmask{seednum}(brainmask));
    
end


for subject = 1:length(subjects)
    subjid = subjects{subject};
    
    outputfilename = ['/data/cn4/evan/RestingState/ConnectivityTimecourses' subjid '.txt'];

delete(outputfilename);
fid = fopen(outputfilename,'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Seed','Network','Window','Value'); %write the output file header
fclose(fid);
dlmwrite(outputfilename,' ','-append');
    
    tmask = dlmread(['/data/cn4/evan/RestingState/28subjects_FDp2_00_DVinit_20_00/tmask_betterDV/' subjid '.txt']);
    tmask = logical(tmask);
    
    data = read_4dfpimg_HCP(['/data/cn4/evan/RestingState/28subjects_FDp2_00_DVinit_20_00/reprocessed_data/' subjid '/' subjid '_333_zmdt_resid_bpss_zmdt_g7.4dfp.img']);
    
    data = data(brainmask,tmask);
    
    
    for seednum = 1:length(seeds)
        
        seed_timecourse(seednum,:) = mean(data(seedmask{seednum},:),1);
    end
        
        
        for windownum = 1:(size(data,2)-windowsize)
            
            string{windownum} = ['Subject ' subjid  ': Window number ' num2str(windownum) ' of ' num2str((size(data,2)-windowsize))];
            if windownum==1; fprintf('%s',string{windownum}); else fprintf([repmat('\b',1,length(string{windownum-1})) '%s'],string{windownum}); end
            
            windowvals = windownum:(windownum+windowsize-1);
            
            wholebraincorrelations = corr(data(:,windowvals)',seed_timecourse(:,windowvals)');
            
            for networknum = 4:numnetworks
                
                for seednum = 1:length(seeds)
                
                avgnetworkcorrelation = mean(wholebraincorrelations(intersect(find(networks==networknum),find(seedmask{seednum}==0)),seednum));
                
                texttowrite = [subjects{subject},'   ',seeds{seednum},'   ',num2str(networknum),'   ',num2str(windownum),'   ',num2str(avgnetworkcorrelation)];  %save the data as a string to be written to the output
                
                dlmwrite(outputfilename,texttowrite,'-append','delimiter','');%write the data to the output file
                
                end
                
            end
        end
        disp('')
        
        clear tmask data seed_timecourse string
    
end

