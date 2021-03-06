%function  surface_parcellation_indsubs_nonstationarity(cohortfile,outputdir,tmasktype,smooth,hems,edges)
%Generate gradient-based parcellation on surface registered subject data
% This script requires a parcellation cohortfile file. The cohortfile includes
% a list of the subject names, the functional volume to be parcellated including
% its path, and the directory where the subject's surface is (as generated by
% Freesurfer and registered to fs_LR space), formatted as:
%
% subjectname funcpath/funcvol surfdir
%
% e.g.
%
% vc33416 /data/cn4/laumannt/vc33416_rest1.4dfp.img /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc33416/7112b_fs_LR
%
% With 'tmasktype', users can enter a tmasklist.txt as created by
% fcimage_analysis to mask out timepoints in their data. Enter 'none' if no
% temporal mask is needed.
% If 'smooth' is set to  'yes' data will be smoothed along the cortical
% surface after surface projection. Enter 'no' if no surface smoothing is
% desired.
% With 'hems', users can specify that either the 'LEFT', 'RIGHT', or 'BOTH'
% hemispheres be parcellated.
% If 'edges' is set to 'yes' a final non-maxima suppression step will be
% applied to the gradient maps and an edge frequency map will be output in
% addition to an average gradient map. Enter 'no' if edge detection is not
% desired.
% This version of surface parcellation will always average correlation data
% from all subjects, generating a single gradient and/or edge map for the
% group.
%Run on tsunami
%matlab12
%addpath /data/cn4/laumannt/assignment_problem_v2/Scripts/
%addpath /data/cn4/laumannt/surface_gradient_code/
%outputdir - full path

% TOL 01/25/13


cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list';
%'/data/cn4/evan/RestingState/FC_Mapping_120/NonstationarityParcellation/LFRS_corrfile2.txt';
surfdatadir = '/data/cn4/evan/RestingState/FC_Mapping_120/NonstationarityParcellation/';
outputdir = '/data/cn4/evan/RestingState/FC_Mapping_120/NonstationarityParcellation/';
tmasktype = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
%'/data/cn4/evan/RestingState/FC_Mapping_120/NonstationarityParcellation/LFRS_tmasklist2.txt';
smooth = 'yes';
hems = 'LEFT';
edges = 'yes';




divisions = 70;
windowsize = 30;
windowoverlap = 0;
smoothnum = 2.55;
mask4dfp='/data/cn4/laumannt/Standard/glm_atlas_mask_333.4dfp.img';
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};

% [ign subnum] = system(['cat ' cohortfile ' | wc -l']);
% subnum = str2num(subnum);
[tmasksubjects tmaskfiles]=textread(tmasktype,'%s%s');
subnum = length(tmasksubjects);

% Read in subject names, functional volume locations, and surface directory
for s = 1:subnum
    [ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
    [ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
    [ign funcvols{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1}''']);
    [ign surfdirs{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' cohortfile]);
end

% Check masktype and the existence of custom masks
% switch tmasktype
%     case 'none'
%     otherwise
%         if ~isequal(tmasksubjects,subjects)
%             error('masklist subjects do not match cohortfile subjects');
%         end
%         for i=1:numel(tmaskfiles)
%             %fprintf('\t%d\t%s\n',i,subjects{i,1});
%             if ~exist(tmaskfiles{i,1})
%                 error('\t%s is not found\n',tmaskfiles{i,1});
%             end
%         end
% end

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

mask = read_4dfpimg(mask4dfp);
        

for hem = h
    
    if numel(subjects)>1
        avgcrosscorr = zeros(32492,65549);
        %avgcrosscorr_surf = zeros(32492,32492);
    end
    maskout = zeros(1,65549);
    
    for s = 100%:length(subjects)
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        subject = strcat(subjects{s});
        funcpath = strcat(funcpaths{s});
        funcvol = strcat(funcvols{s});
        surfdir = strcat(surfdirs{s});
        
        % Move funcvol to outputdir and convert to nifti
        evalc(['!cp ' funcpath '/' funcvol '.4dfp.* ' outputdir]);
        cd(outputdir)
%         evalc(['!nifti_4dfp -n ' outputdir '/' funcvol '.4dfp.img ' outputdir '/' funcvol '.nii']);
%         
%         % Sample timecourse to surface
%                 disp('Projecting timecourses to surface')
        BOLDname = [subject '_BOLD_' HEMS{hem}];
        
        
%         if strcmp(subject(end-1),'_')
%             midsurf = [surfdir '/Native/' subject(1:end-2) '.' HEMS{hem} '.midthickness.native.surf.gii'];
%             whitesurf = [surfdir '/Native/' subject(1:end-2) '.' HEMS{hem} '.white.native.surf.gii'];
%             pialsurf = [surfdir '/Native/' subject(1:end-2) '.' HEMS{hem} '.pial.native.surf.gii'];
%         else
%             midsurf = [surfdir '/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
%             whitesurf = [surfdir '/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
%             pialsurf = [surfdir '/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
%         end
%         
%         system([workbenchdir '/wb_command -volume-to-surface-mapping ' outputdir '/' funcvol '.nii ' midsurf ' ' BOLDname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf ' -voxel-subdiv 5']);
%         
%         % Smooth BOLD data on surface
        switch smooth
            case 'no'
            case 'yes'
%                 disp('Smoothing data on surface')
%                 system([workbenchdir '/wb_command -metric-smoothing ' midsurf ' ' BOLDname '.func.gii ' num2str(smoothnum) ' ' BOLDname '_smooth' num2str(smoothnum) '.func.gii'])
%                 system(['rm ' BOLDname '.func.gii'])
                BOLDname = [BOLDname '_smooth' num2str(smoothnum)];
        end
        
        % Deform from native to 32K fs_LR surface
%         disp('Deform timecourse to 32k fs_LR')
%         cd([surfdir '/fsaverage_LR32k/'])
%         system(['caret_command64 -deformation-map-apply native232k_fs_LR.' HEMS{hem} '.deform_map METRIC_AVERAGE_TILE ' outputdir '/' BOLDname '.func.gii ' outputdir '/' BOLDname '_32k_fsLR.func.gii']);
%         
%         % Convert surface timecourses to be read by matlab
         cd(surfdatadir)
%         system(['rm ' BOLDname '.func.gii'])
        BOLDname = [BOLDname '_32k_fsLR'];
%         system(['caret_command64 -file-convert -format-convert ASCII ' BOLDname '.func.gii'])
        evalc(['!rm ' BOLDname '_noHEAD.func.gii']);
        evalc(['!awk ''NF > 25'' ' BOLDname '.func.gii > ' BOLDname '_noHEAD.func.gii']);
        
        % Load surface and volume data
        surf_BOLD = load([BOLDname '_noHEAD.func.gii']);
        surf_BOLD(:,1) = [];
        
        cd(outputdir)
        
        BOLD = read_4dfpimg(fullfile(outputdir, [funcvol '.4dfp.img']));
        vol_BOLD = BOLD(logical(mask),:)';
        
        % Tmask data
        switch tmasktype
            case 'none'
            otherwise
                tmask = load(tmaskfiles{s});
                surf_BOLD = surf_BOLD(:,logical(tmask));
                vol_BOLD = vol_BOLD(logical(tmask),:);
        end
        
        
        
        numwindows = floor((size(surf_BOLD,2)-windowsize)/(windowsize-windowoverlap))+1;
        %crosscorr = zeros(size(surf_BOLD,1),size(vol_BOLD,2)*numwindows,'single');
        
        %windowsperdivision = 5;
        %divisions = numwindows * windowsperdivision;
        
        nodesperdivision = ceil(size(surf_BOLD,1) / divisions);
        
        %disp('Calculating correlation between surface nodes and volume')
        for divisionnum = 1:divisions
            
            if divisionnum==divisions
                indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : size(surf_BOLD,1);
            else
                indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : nodesperdivision*divisionnum;
            end
            
        end
%             
%             crosscorr = zeros(length(indices{divisionnum}),size(vol_BOLD,2)*numwindows,'single');
%             
%             for windownum = 1:numwindows
%                 % Calculate correlation between surface and volume data
%                 
%                 string{windownum} = ['  Division ' num2str(divisionnum) ' of ' num2str(divisions) ', window ' num2str(windownum) ' of ' num2str(numwindows)];
%                 if windownum==1; fprintf('%s',string{windownum}); else fprintf([repmat('\b',1,length(string{windownum-1})) '%s'],string{windownum}); end
%                 
%                 windowindices = [(windowsize-windowoverlap)*(windownum-1)+1 : (windowsize-windowoverlap)*(windownum-1) + windowsize];
%                 
%                 crosscorr(:,size(vol_BOLD,2)*(windownum-1)+1 : size(vol_BOLD,2)*windownum) = FisherTransform(paircorr_mod(single(surf_BOLD(indices{divisionnum},windowindices))',single(vol_BOLD(windowindices,:))));
%                 
%                 if windownum==1 && divisionnum==1
%                     medialmaskcheckmatrix = FisherTransform(paircorr_mod(single(surf_BOLD)',single(vol_BOLD)));
%                     medialmaskcheckmatrix = medialmaskcheckmatrix(:,1);
%                 end
%                 
%             end
%             
%             
%             disp(' ')
%             
%             % Remove voxels from crosscorr that have no BOLD data but are in glm mask
%             nanmat = isnan(crosscorr);
%             nancolsum = sum(nanmat,2)==size(crosscorr,2);
%             nanmat(logical(nancolsum),:) = [];
%             nanrowsum = sum(nanmat,1)>0;
%             %nanrowsum = repmat(nanrowsum,1,numwindows);
%             crosscorr(:,logical(nanrowsum)) = [];
%             
%             save([outputdir '/CrosscorrDivision' num2str(divisionnum) '.mat'],'crosscorr','-v7.3')
%             
%             clear crosscorr
%             
%         end
        
        %crosscorr_surf = FisherTransform(paircorr_mod(single(surf_BOLD)'));
        
        % Keep running correlation map average
        %         if numel(subjects)>1
        %             avgcrosscorr = avgcrosscorr + crosscorr;
        %             %avgcrosscorr_surf = avgcrosscorr_surf + crosscorr_surf;
        %         end
        
        % Keep track of voxels without data by subject
        %maskout = maskout + isnan(crosscorr(1,size(vol_BOLD,2)));
        
        evalc(['!rm ' BOLDname '_noHEAD.func.gii']);
        evalc(['!rm ' funcvol '.*']);
        
        
        
        %         avgcrosscorr = crosscorr;
        %         clear crosscorr
        
        
        
        
        
        %save(gifti(single(medialmask)),'medialmask.func.gii')
        
        %for piecenum = 1:(numwindows*2)
        
        % Calculate correlation of correlation maps
        fprintf('Calculating correlation of correlation maps...')
        disp(' ')
        avgcorrofcorr = zeros(size(surf_BOLD,1),'single');
        
        
        
        divisionnum = 1;
        stringcounter = 1;
        
        for i = 1:divisions
            
            crosscorri = zeros(length(indices{i}),size(vol_BOLD,2)*numwindows,'single');
            
            for windownum = 1:numwindows
                % Calculate correlation between surface and volume data
                
                windowindices = [(windowsize-windowoverlap)*(windownum-1)+1 : (windowsize-windowoverlap)*(windownum-1) + windowsize];
                
                crosscorri(:,size(vol_BOLD,2)*(windownum-1)+1 : size(vol_BOLD,2)*windownum) = FisherTransform(paircorr_mod(single(surf_BOLD(indices{i},windowindices))',single(vol_BOLD(windowindices,:))));
                
                if windownum==1 && divisionnum==1
                    medialmaskcheckmatrix = FisherTransform(paircorr_mod(single(surf_BOLD)',single(vol_BOLD)));
                    medialmaskcheckmatrix = medialmaskcheckmatrix(:,1);
                end
            end
                
            % Remove voxels from crosscorr that have no BOLD data but are in glm mask
            nanmat = isnan(crosscorri);
            nancolsum = sum(nanmat,2)==size(crosscorri,2);
            nanmat(logical(nancolsum),:) = [];
            nanrowsum = sum(nanmat,1)>0;
            %nanrowsum = repmat(nanrowsum,1,numwindows);
            crosscorri(:,logical(nanrowsum)) = [];
                
            
            
            
            for j = 1:divisions
                
                if i >= j
                    
                    
                    
                    crosscorrj = zeros(length(indices{j}),size(vol_BOLD,2)*numwindows,'single');
                    
                    for windownum = 1:numwindows
                        % Calculate correlation between surface and volume data
                        
                        string{stringcounter} = ['Calculating correlations for division number ' num2str(divisionnum) ' of ' num2str((divisions^2 + divisions)/2) ': window ' num2str(windownum) ' of ' num2str(numwindows)];
                        if stringcounter==1; fprintf('%s',string{stringcounter}); else fprintf([repmat('\b',1,length(string{stringcounter-1})) '%s'],string{stringcounter}); end
                        stringcounter = stringcounter + 1;
                    
                        
                        windowindices = [(windowsize-windowoverlap)*(windownum-1)+1 : (windowsize-windowoverlap)*(windownum-1) + windowsize];
                        
                        crosscorrj(:,size(vol_BOLD,2)*(windownum-1)+1 : size(vol_BOLD,2)*windownum) = FisherTransform(paircorr_mod(single(surf_BOLD(indices{j},windowindices))',single(vol_BOLD(windowindices,:))));
                        

                    end
                    
                    % Remove voxels from crosscorr that have no BOLD data but are in glm mask
                    nanmat = isnan(crosscorrj);
                    nancolsum = sum(nanmat,2)==size(crosscorrj,2);
                    nanmat(logical(nancolsum),:) = [];
                    nanrowsum = sum(nanmat,1)>0;
                    %nanrowsum = repmat(nanrowsum,1,numwindows);
                    crosscorrj(:,logical(nanrowsum)) = [];
                    
                    avgcorrofcorr(indices{i},indices{j}) = paircorr_mod(crosscorri',crosscorrj');
                    
                    divisionnum = divisionnum+1;
                    
                end
                
            end
        end
        disp(' ')
        
        clear crosscorri crosscorrj
        
        % Create mask based on nodes that don't get projected to
        medialind = isnan(medialmaskcheckmatrix);
        medialmask = ones(length(medialmaskcheckmatrix),1);
        medialmask(medialind) = 0;
        
        avgcorrofcorr = tril(avgcorrofcorr,-1) + tril(avgcorrofcorr)';
        
        system(['rm ' outputdir '/CrosscorrDivision*.mat']);
        
        disp('done.')
        avgcorrofcorr(isnan(avgcorrofcorr)) = 0;
        avgcorrofcorr(medialind,:) = [];
        
        % Save metrics
        cd(outputdir)
        fprintf('Saving corr of corr maps...')
        save(gifti(avgcorrofcorr(1:floor(size(avgcorrofcorr,1)/2),:)),['avg_corrofcorr_' HEMS{hem} '_piece1.func.gii'],'ExternalFileBinary');
        save(gifti(avgcorrofcorr(floor(size(avgcorrofcorr,1)/2)+1:end,:)),['avg_corrofcorr_' HEMS{hem} '_piece2.func.gii'],'ExternalFileBinary');
        disp('done.')
        
        clear avgcorrofcorr
        
        % Calculate gradient
        atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser';
        specfile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.32k_fs_LR.c5.spec'];
        midsurf_32k = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        coordfile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.coord.gii'];
        topofile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.32k_fs_LR.topo.gii'];
        
        fprintf('Calculating Gradients...')
        filestem = [subject '_avgcorrofcorr_smooth' num2str(smoothnum) '_allgrad_' HEMS{hem} '_nonstationarityORIGDATA'];
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
%         grad_metrics_piece1 = load([outputdir '/' filestem '_piece1_noHEAD.func.gii']);
%         grad_metrics_piece1(:,1) = [];
%         grad_metrics_piece2 = load([outputdir '/' filestem '_piece2_noHEAD.func.gii']);
%         grad_metrics_piece2(:,1) = [];
%         grad_metrics = [grad_metrics_piece1 grad_metrics_piece2];
%         
%         % Calculate and save average gradient map
%         avg_metric = mean(grad_metrics,2);
%         save(gifti(single(avg_metric)),[filestem '_avg.func.gii']);
%         
%         % Calculate edge maps
%         switch edges
%             case 'no'
%             case 'yes'
%                 surface_edges_all_test_faster(grad_metrics,specfile,outputdir,[filestem '_edge_avg'],1);
%         end
        
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
        
        clear avg_metric grad_metrics grad_metrics_piece1 grad_metrics_piece2
        
    end
end
