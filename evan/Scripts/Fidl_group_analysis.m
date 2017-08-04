warning off

subjects = {'vc34096' 'vc34116' 'vc34125' 'vc34126' 'vc34128' 'vc34140' 'vc34141' 'vc34198' 'vc34199' 'vc34200' 'vc34201' 'vc34220' 'vc34250' 'vc34251' 'vc34252' 'vc34253' 'vc34268' 'vc34306' 'vc34307' 'vc34308' 'vc34330' 'vc34331' 'vc34401' 'vc34402' 'vc34403' 'vc34404' 'vc34405' 'vc34408'};

timepoints = 7;

voxelspace = '333';

outputname = 'METaskxTime';

linuxdir = '/home/usr/fidl/fidl_code/fidl_2.64/bin_linux';

factornames = {'METask'};

levelnames = {'Glass','MRT','Semantic'};

ANOVAsetup = {[4 5 6 7],[13 14 15 16 17 18],[24 25]};

basedir = '/data/cn3/joe/ResourceDataLimited/';

GLMsuffix = '_Scrubbed_f06_NameShortened';

T4suffix = '_anat_ave_to_TRIO_KY_NDC_t4';

Eventfilesuffix = '_EventFile_NameShortened';

outputdir = '/data/cn4/evan/Occipitotemporal/Joe3Tasks/';

Scratchdir = 'SCRATCH/';

RunSingleSubs = 1;

mkdir(outputdir);
mkdir([outputdir Scratchdir]);
if RunSingleSubs
    mkdir([outputdir '/SingleSubjects/']);
end

driverfile{1} = ['subject    '];
    for factor = 1:length(factornames)
        driverfile{1} = [driverfile{1}  factornames{factor}];
    end
    
    if numel(ANOVAsetup) > 1
        driverfile{1} = [driverfile{1} '    time'];
    end
    
    driverfile{1} = [driverfile{1}  '    *.4dfp.img'];

    driverfilename = [outputdir outputname '_driver.dat'];
    
delete(driverfilename);
fid = fopen(driverfilename,'at'); %open the output file for writing
fprintf(fid,'%s\n',driverfile{1}); %write the output file header
fclose(fid);
 
    
anovarunstr = ['!' linuxdir '/fidl_anova_new -driver "' driverfilename '" -voxel_threshold 0.01 -uncompress /home/usr/fidl/lib/glm_atlas_mask_' voxelspace '.4dfp.img -output Z_uncorrected Z_monte_carlo  -glmfiles '];



for subject = 1:length(subjects)
    
    clear conditionindex
    
    subname = subjects{subject};
    
    disp(['Subject ' subname])
    
    
    rawonsets = dlmread([basedir subname '/' subname Eventfilesuffix '.fidl'],'\t',1,0);
    onsets = rawonsets(:,1:3);
    rawstuff = textread([basedir subname '/' subname Eventfilesuffix '.fidl'],'%s\t');
    TR = str2num(rawstuff{1});
    stilltext = 1;
    entrynum = 2;
    while stilltext==1
        if isempty(str2num(rawstuff{entrynum}))
            conditionnames{entrynum-1} = rawstuff{entrynum};
            entrynum = entrynum+1;
        else
            stilltext = 0;
        end
    end
    
    
    for condnum = 1:length(conditionnames)
        
        thisonsetindices = find(onsets(:,2) == condnum-1);
        
        if condnum == 1
            conditionindex(condnum) = 1;
            
            if onsets(thisonsetindices(1),3) > 2
                    
                    lastduration = 1;
                else
                    lastduration = timepoints;
                end
        else
            
            
            
            
            if ~isempty(thisonsetindices)
                
                conditionindex(condnum) = conditionindex(condnum-1) + lastduration;
            
                
                if onsets(thisonsetindices(1),3) > 2
                    
                    lastduration = 1;
                else
                    lastduration = timepoints;
                end
            
            else
                
                conditionindex(condnum) = conditionindex(condnum-1);
                
            end
            
             
        end
    end
        
        
    subjzstatrunstr = ['!' linuxdir '/fidl_zstat -glm_file ' basedir subname '/' subname GLMsuffix '.glm -tc '];
    singlesubjrunstr = ['!' linuxdir '/fidl_avg_zstat -glm_files ' basedir subname '/' subname GLMsuffix '.glm -tc '];
    
     for cell = 1:numel(ANOVAsetup)
         string = '[';
         if isempty(find(size(ANOVAsetup)==1))
            nfactors = ndims(ANOVAsetup);
         else
             nfactors = 1;
         end
         
         for factor = 1:nfactors;
             string = [string 'dim{' num2str(factor) '} '];
         end
         string = [string '] = ind2sub(size(ANOVAsetup),cell)'];
         evalc(string);
         
         
         singlesubjectoutfilename{cell} = ['"' subname '_avgtc_' levelnames{cell}];

            
            for time = 1:timepoints
                
                singlesubjectoutfilename{cell} = [singlesubjectoutfilename{cell} '_' num2str(time)];
                
                timecourse{subject}{time}{cell} = [];
                
                driverfile{end+1} = num2str(subject); 
                
                if numel(ANOVAsetup)>1
                    
                    string = [];
                    
                for factor = 1:length(dim)
                    
                    string = [string num2str(dim{factor}) ','];
                    
                end
                
                string = string(1:end-1);
                
                    evalc(['driverfile{end} = [driverfile{end} ''    '' levelnames{' string '}]']);
                
                end
                
                driverfile{end} = [driverfile{end}  '    ' num2str(time)];

                driverfile{end} = [driverfile{end} '    ' Scratchdir subname GLMsuffix '_'];
                for condnum = 1:length(ANOVAsetup{cell})
                    
                    driverfile{end} = [driverfile{end} conditionnames{ANOVAsetup{cell}(condnum)} '+'];
                    
                    timecourse{subject}{time}{cell} = [timecourse{subject}{time}{cell} num2str(conditionindex(ANOVAsetup{cell}(condnum)) +time -1) '+'];
                    
                end
                timecourse{subject}{time}{cell} = [timecourse{subject}{time}{cell}(1:end-1) ' '];
                
                subjzstatrunstr = [subjzstatrunstr timecourse{subject}{time}{cell} ' '];
                singlesubjrunstr = [singlesubjrunstr timecourse{subject}{time}{cell} ' '];
                
                
                driverfile{end} = [driverfile{end}(1:end-1) '_' num2str(time) '.4dfp.img'];
                
                dlmwrite(driverfilename,driverfile{end},'-append','delimiter','');
                
            end
            
            
            singlesubjectoutfilename{cell} = [singlesubjectoutfilename{cell} '"'];
                
            
     end
     
     
     
     subjzstatrunstr = [subjzstatrunstr '-xform_file ' basedir subname '/atlas/' subname T4suffix ' -gauss_smoth 2 -compress /home/usr/fidl/lib/glm_atlas_mask_' voxelspace '.4dfp.img -atlas ' voxelspace ' -scratchdir ' outputdir Scratchdir];
     singlesubjrunstr = [singlesubjrunstr '-xform_files ' basedir subname '/atlas/' subname T4suffix ' -mask /home/usr/fidl/lib/glm_atlas_mask_' voxelspace '.4dfp.img -atlas ' voxelspace ' -directory ' outputdir '/SingleSubjects/  -frames ' num2str(ones(1,numel(ANOVAsetup))*timepoints) ' -glmpersub 1 -tc_names '];
     
     for cell = 1:numel(ANOVAsetup)
         singlesubjrunstr = [singlesubjrunstr singlesubjectoutfilename{cell} ' '];
     end
     
     sublog{subject} = evalc(subjzstatrunstr);
     
     if RunSingleSubs
         singlesublog{subject} = evalc(singlesubjrunstr);
     end
     
     
     anovarunstr = [anovarunstr basedir subname '/' subname GLMsuffix '.glm '];
    
end

anovarunstr = [anovarunstr '-tc '];

for cell = 1:numel(ANOVAsetup)
for time = 1:timepoints
    for subject = 1:length(subjects)
        anovarunstr = [anovarunstr timecourse{subject}{time}{cell}(1:end-1) ','];
    end
    anovarunstr = [anovarunstr(1:end-1) ' '];
end
end

cd(outputdir)

disp('Group ANOVA')

anovarunstr = [anovarunstr '-Nimage_name "' outputdir outputname '_Nimage.4dfp.img" -threshold_extent "3.75 18" -pval .05 -scratchdir ' outputdir Scratchdir ' -GIGAdesign -glmpersub ' num2str(ones(1,length(subjects))) ' -clean_up'];
anovalog = evalc(anovarunstr);
     




    