%Generate gradient-based parcellation on surface registered subject beta
%images


% TOL 01/25/13

outputdir = '/data/cn4/evan/Task_parcellation/Stevedata_subregister/';
datadir = '/data/cn4/evan/Task_parcellation/Manytasks_subregister/';
smooth = 'yes';
hems = 'BOTH';
edges = 'yes';


divisions = 6;


lists = {'/data/cn4/evan/Task_parcellation/SteveSubjects/SteveFakeList_newdata.glm_list'};
%     '/data/cn4/evan/Occipitotemporal/VisualSwitch/JessTaskSwitchCueing_ANC_CsepT_SCRUBfbf07_all30_9tps_020312.glm_list',...
%     '/data/cn4/evan/Occipitotemporal/VisualAttention/fbf05_GLMs.txt',...
%     '/data/cn4/evan/28Subjects_RePreprocessed_Scrubbed_f06_NameShortened_wRestData.list',...
 


directories = {'/data/cn4/evan/Task_parcellation/SteveSubjects/OrigData/'};
%     '/data/cn4/evan/Occipitotemporal/VisualSwitch/CsepT_Trials/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/VisualAttention/SingleSubjects/',...
%     '/data/cn4/evan/Occipitotemporal/Joe3Tasks_IndConditions/SingleSubjects/',...
    


smoothnum = 2.55;
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};


%Read the brain mask
mask4dfp='/data/cn4/laumannt/Standard/glm_atlas_mask_333.4dfp.img';
volmask = read_4dfpimg(mask4dfp);


switch hems
    case 'LEFT'
        h = 1;
    case 'RIGHT'
        h = 2;
    case 'BOTH'
        h = [1:2];
end

% Make output folder
evalc(['!mkdir ' outputdir]);


for hem = h
    
    medialmask = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
    medialmaskdata = gifti(medialmask);
    medialind = find(medialmaskdata.cdata);
    
    %figure out the size of each division
    nodesperdivision = ceil(32492 / divisions);
    
    warning off
    
    disp('Calculating correlation between surface nodes and volume')
    
    %Loop through divisions
    for divisionnum = 1:divisions
        numallsubjects = 0;
        %get the vertex indices for this division
        if divisionnum==divisions
            indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : 32492;
        else
            indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : nodesperdivision*divisionnum;
        end
        subcounter = 0;
        for listnum = 1:length(lists)
            list = lists{listnum};
            %disp(list)
            directory = directories{listnum};
            %timepoints = listtimepoints(listnum);
            %extranamestuff = extranamestuffstrings{listnum};
            
            [front entries] = textread(list,'%s%s','delimiter',':');
            
            numsubjects = str2num(entries{2});
            glmfiles = entries(3:2+numsubjects);
            %T4files = entries(4+numsubjects:end);
            for i = 1:length(glmfiles)
                if strcmp(glmfiles{i},'')
                    glmfiles{i} = front{2+i};
                end
                %     if strcmp(T4files{i},'')
                %         T4files{i} = front{3+numsubjects+i};
                %     end
            end
            
            numallsubjects = numallsubjects + length(glmfiles);
            
            %Loop through subjects
            for s = 1:length(glmfiles)
                
                try
                    subnameindex = strfind(glmfiles{s},'/vc');
                    slashindex = min([strfind(glmfiles{s}(subnameindex(1)+1:end),'/') strfind(glmfiles{s}(subnameindex(1)+1:end),'_') strfind(glmfiles{s}(subnameindex(1)+1:end),'.')]);
                    subname = glmfiles{s}(subnameindex+1:subnameindex+slashindex-1);
                    conditionfiles = dir([directory subname '*.4dfp.img']);
                catch
                    subname = [glmfiles{s}];
                    conditionfiles = dir([directory '_' subname '*.4dfp.img']);
                end
                subcounter = subcounter+1;
                string{subcounter} = ['    Division ' num2str(divisionnum) ' of ' num2str(divisions) ', subject #' num2str(subcounter) ': ' subname];
                if s==1; fprintf('%s',string{subcounter}); else fprintf([repmat('\b',1,length(string{subcounter-1})) '%s'],string{subcounter}); end
                
                voldata = [datadir '/' subname '_Mergefile.4dfp.img'];
                surfdata = [voldata(1:end-9) '_' HEMS{hem} '_32k_fsLR.func.gii'];
                if divisionnum == 1
                    
                    if ~exist(voldata)
                        
                        fslmergestr = ['fslmerge -t ' voldata(1:end-9) '.nii.gz'];
                        
                        for conditionnum = 1:length(conditionfiles)
                            fslmergestr = [fslmergestr ' ' directory conditionfiles(conditionnum).name(1:end-9) '.nii'];
                            
                            if ~exist([directory conditionfiles(conditionnum).name(1:end-9) '.nii'])
                                evalc(['!nifti_4dfp -n ' directory conditionfiles(conditionnum).name ' ' directory conditionfiles(conditionnum).name(1:end-9) '.nii']);
                            end
                            
                        end
                        
                        evalc(['!' fslmergestr]);
                        
                        clear fslmergestr
                        
                        gunzip([voldata(1:end-9) '.nii.gz']);
                        evalc(['!nifti_4dfp -4 ' voldata(1:end-9) '.nii ' voldata ]);
                        
                    end
                    
                    if ~exist(surfdata)
                        
                        surfdir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subname '/7112b_fs_LR/'];
                        
                        midsurf = [surfdir '/Native/' subname '.' HEMS{hem} '.midthickness.native.surf.gii'];
                        whitesurf = [surfdir '/Native/' subname '.' HEMS{hem} '.white.native.surf.gii'];
                        pialsurf = [surfdir '/Native/' subname '.' HEMS{hem} '.pial.native.surf.gii'];
                        
                        fileexists = 0;
                        while fileexists == 0
                            fileexists = exist(midsurf);
                            pause(10)
                        end
                        
                        evalc(['!' workbenchdir '/wb_command -volume-to-surface-mapping '  voldata(1:end-9) '.nii ' midsurf ' ' voldata(1:end-9) '_' HEMS{hem} '_native.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -voxel-subdiv 5']);
                        evalc(['!' workbenchdir '/wb_command -metric-smoothing ' midsurf ' ' voldata(1:end-9) '_' HEMS{hem} '_native.func.gii ' num2str(smoothnum) ' ' voldata(1:end-9) '_' HEMS{hem} '_native_smooth.func.gii']);
                        
                        cd([surfdir '/fsaverage_LR32k/'])
                        evalc(['!caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE '  voldata(1:end-9) '_' HEMS{hem} '_native_smooth.func.gii '  surfdata]);
                        cd(outputdir)
                        
                        evalc(['!caret_command64 -file-convert -format-convert XML ' surfdata]);
                        
                        delete([outputdir '/' subname '_Mergefile.nii.gz'])
                        delete([outputdir '/' subname '_Mergefile_' HEMS{hem} '_native.func.gii'])
                        delete([outputdir '/' subname '_Mergefile_' HEMS{hem} '_native_smooth.func.gii'])
                    end
                end
                
                %Read in the surface data and get the betas to be parcellated
                %surf_BOLD = load([surfdata(1:end-9) '_noHEAD.func.gii']);
                %surf_BOLD(:,1) = [];
                surf_BOLD = gifti(surfdata);
                surf_BOLD = surf_BOLD.cdata;
                
                %Read in the volume data and get the betas to be parcellated
                %voldata4dfp = [voldata(1:end-7) '.4dfp.img'];
                BOLD = read_4dfpimg(voldata);
                vol_BOLD = BOLD(logical(volmask),:)';
                vol_BOLD(find(vol_BOLD==0)) = NaN;
                
                %Set up the avgcrosscorr variable (avg across subjects) for this division
                if s==1 && listnum==1
                    avgcrosscorr = zeros(length(indices{divisionnum}), size(vol_BOLD,2));
                end
                
                %Calculate the surface-volume correlation
                crosscorr = FisherTransform(paircorr_mod(single(surf_BOLD(indices{divisionnum},:))',single(vol_BOLD)));
                
                %add that correlation to the avgcrosscorr variable
                avgcrosscorr = avgcrosscorr + crosscorr;
                
                clear crosscorr
                
            end
        end
        
        %Track progress
        string{end+1} = ['    Division ' num2str(divisionnum) ' of ' num2str(divisions) ' complete'];
        fprintf([repmat('\b',1,length(string{subcounter})) '%s'],string{end});
        clear string
        disp(' ')
        
        %Divide the avgcrosscorr by the number of subjects
        avgcrosscorr = avgcrosscorr ./ numallsubjects;
        
        % Remove voxels from crosscorr that have no BOLD data but are in glm mask
        nanmat = isnan(avgcrosscorr);
        nancolsum = sum(nanmat,2)==length(avgcrosscorr(1,:));
        nanmat(logical(nancolsum),:) = [];
        nanrowsum = sum(nanmat,1)>0;
        avgcrosscorr(:,logical(nanrowsum)) = [];
        
        %Save the avgcrosscorr for this division
        save([outputdir '/CrosscorrDivision_' num2str(divisionnum) '.mat'],'avgcrosscorr','-v7.3')
        
        clear avgcrosscorr string
        
    end
    
    %Delete all the converted subject data (to save disk space)
    %         delete([outputdir '/*_Mergefile.4dfp.img'])
    %         delete([outputdir '/*_' HEMS{hem} '*.func.gii'])
    
    disp('done.')
    
    
    fprintf('Calculating correlation of correlation maps...')
    disp(' ')
    avgcorrofcorr = zeros(size(surf_BOLD,1),'single');
    divisionnum = 0;
    
    %Loop through division pairs
    for i = 1:divisions
        
        %Load data from one division
        load([outputdir '/CrosscorrDivision_' num2str(i) '.mat']);
        crosscorri = avgcrosscorr;
        
        for j = 1:divisions
            
            if i >= j
                
                %Load data from another division
                load([outputdir '/CrosscorrDivision_' num2str(j) '.mat']);
                crosscorrj = avgcrosscorr;
                
                %Track progress
                divisionnum = divisionnum+1;
                string{divisionnum} = ['Calculating correlations for division number ' num2str(divisionnum) ' of ' num2str((divisions^2 + divisions)/2)];
                if divisionnum==1; fprintf('%s',string{divisionnum}); else fprintf([repmat('\b',1,length(string{divisionnum-1})) '%s'],string{divisionnum}); end
                
                %Calculate the corr of corr between these two divisions
                avgcorrofcorr(indices{i},indices{j}) = paircorr_mod(crosscorri',crosscorrj');
                
            end
            
        end
    end
    disp(' ')
    
    clear avgcrosscorr crosscorri crosscorrj
    
    %We've only calculated the lower triangular half of this matrix, so create the reflection and add it to variable to complete the matrix
    avgcorrofcorr = tril(avgcorrofcorr,-1) + tril(avgcorrofcorr)';
    
    %Delete all the temporary division variables
    %system(['rm ' outputdir '/CrosscorrDivision_*.mat']);
    
    %         avgcorrofcorr = paircorr_mod(single(avgcrosscorr)');
    %         clear avgcrosscorr
    cd(outputdir)
    
    %         avgcorrofcorrStevedata = avgcorrofcorr;
    %         save([outputdir '/avgcorrofcorrStevedata.mat'],'avgcorrofcorrStevedata','-v7.3')
    %         clear all
    %         surface_taskparcellation_bysubs
    
    
    avgcorrofcorr(isnan(avgcorrofcorr)) = 0;
    avgcorrofcorr(medialind,:) = [];
    
    % Save metrics
    cd(outputdir)
    fprintf('Saving corr of corr maps...')
    save(gifti(avgcorrofcorr(1:floor(size(avgcorrofcorr,1)/2),:)),['avg_corrofcorr_' HEMS{hem} '_piece1.func.gii'],'ExternalFileBinary');
    save(gifti(avgcorrofcorr(floor(size(avgcorrofcorr,1)/2)+1:end,:)),['avg_corrofcorr_' HEMS{hem} '_piece2.func.gii'],'ExternalFileBinary');
    disp('done.')
    
    % Calculate gradient
    atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser';
    specfile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.32k_fs_LR.c5.spec'];
    midsurf_32k = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    coordfile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.coord.gii'];
    topofile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.32k_fs_LR.topo.gii'];
    
    
    
    fprintf('Calculating Gradients...')
    filestem = ['avgcorrofcorr_smooth' num2str(smoothnum) '_allgrad_' HEMS{hem} ];
    system(['caret_command64 -metric-gradient-all ' coordfile ' ' topofile ' avg_corrofcorr_' HEMS{hem} '_piece1.func.gii ' filestem '_piece1.func.gii FALSE ' num2str(smoothnum) ' TRUE']);
    system(['caret_command64 -metric-gradient-all ' coordfile ' ' topofile ' avg_corrofcorr_' HEMS{hem} '_piece2.func.gii ' filestem '_piece2.func.gii FALSE ' num2str(smoothnum) ' TRUE']);
    
    system(['rm avg_corrofcorr_' HEMS{hem} '_piece*'])
    
    % Reformatting files to load into Matlab
    fprintf('Reformatting gifti files to load into matlab')
    system(['caret_command64 -file-convert -format-convert ASCII ' outputdir '/' filestem '_piece1.func.gii'])
    system(['caret_command64 -file-convert -format-convert ASCII ' outputdir '/' filestem '_piece2.func.gii'])
    system(['rm ' outputdir '/' filestem '_piece1_noHEAD.func.gii'])
    system(['rm ' outputdir '/' filestem '_piece2_noHEAD.func.gii'])
    system(['awk ''NF > 25'' ' outputdir '/' filestem '_piece1.func.gii > ' filestem '_piece1_noHEAD.func.gii'])
    system(['awk ''NF > 25'' ' outputdir '/' filestem '_piece2.func.gii > ' filestem '_piece2_noHEAD.func.gii'])
    
    
    % Extra smoothing step...
    fprintf('Smoothing gradients before calculating edges...')
    system(['caret_command64 -metric-smoothing ' coordfile ' ' topofile ' ' filestem '_piece1.func.gii ' filestem '_smooth' num2str(smoothnum) '_piece1.func.gii GEOGAUSS 1 1 -geo-gauss ' num2str(smoothnum)])
    system(['caret_command64 -metric-smoothing ' coordfile ' ' topofile ' ' filestem '_piece2.func.gii ' filestem '_smooth' num2str(smoothnum) '_piece2.func.gii GEOGAUSS 1 1 -geo-gauss ' num2str(smoothnum)])
    
    % Remove intermediate file
    system(['rm ' filestem '_piece*']);
    system(['rm ' filestem '_piece*']);
    
    % Reformatting files to load into Matlab
    filestem = [filestem '_smooth' num2str(smoothnum)];
    fprintf('Reformatting gifti files to load into matlab')
    system(['caret_command64 -file-convert -format-convert ASCII ' outputdir '/' filestem '_piece1.func.gii'])
    system(['caret_command64 -file-convert -format-convert ASCII ' outputdir '/' filestem '_piece2.func.gii'])
    system(['rm ' outputdir '/' filestem '_piece1_noHEAD.func.gii'])
    system(['rm ' outputdir '/' filestem '_piece2_noHEAD.func.gii'])
    system(['awk ''NF > 25'' ' outputdir '/' filestem '_piece1.func.gii > ' filestem '_piece1_noHEAD.func.gii'])
    system(['awk ''NF > 25'' ' outputdir '/' filestem '_piece2.func.gii > ' filestem '_piece2_noHEAD.func.gii'])
    grad_metrics_piece1 = load([outputdir '/' filestem '_piece1_noHEAD.func.gii']);
    grad_metrics_piece1(:,1) = [];
    grad_metrics_piece2 = load([outputdir '/' filestem '_piece2_noHEAD.func.gii']);
    grad_metrics_piece2(:,1) = [];
    grad_metrics = [grad_metrics_piece1 grad_metrics_piece2];
    
    % Clean up files
    system(['rm ' filestem '_piece*']);
    system(['rm ' filestem '_piece*']);
    
    % Calculate and save average gradient map
    avg_metric = mean(grad_metrics,2);
    save(gifti(avg_metric),[filestem '_avg.func.gii']);
    
    % Calculate edge maps
    switch edges
        case 'no'
        case 'yes'
            surface_edges_all_test_faster(grad_metrics,specfile,outputdir,[filestem '_edge_avg'],1);
    end
    
    
    
    
end
%end
