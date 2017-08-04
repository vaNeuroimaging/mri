%Generate gradient-based parcellation on surface registered subject beta
%images


% TOL 01/25/13

outputdir = '/data/cn4/evan/Task_parcellation/SteveSubjects/Parcellation_newdata_bysubs/';
datadir = '/data/cn4/evan/Task_parcellation/SteveSubjects/';
smooth = 'no';
hems = 'LEFT';
edges = 'yes';
medialmask = '/data/cn4/evan/ROIs/medialwall.func.gii';
medialmaskdata = gifti(medialmask);
medialind = find(medialmaskdata.cdata);


divisions = 6;


subjectlist = '/data/cn4/evan/Task_parcellation/SteveSubjects/Subjectlist2.txt';

subjects = textread(subjectlist,'%s');

    
    smoothnum = 2.55;
    workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
    HEMS = {'L';'R'};
    hemname = {'LEFT';'RIGHT'};
    
    hems = 'LEFT';
    
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
        
        cd(datadir)
        
        %figure out the size of each division
        nodesperdivision = ceil(32492 / divisions);
        
        warning off
        
        disp('Calculating correlation between surface nodes and volume')
        
        %Loop through divisions
        for divisionnum = 1:divisions
            
            %get the vertex indices for this division
            if divisionnum==divisions
                indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : 32492;
            else
                indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : nodesperdivision*divisionnum;
            end
            
            %Lopp thhrough subjects
            for s = 1:length(subjects)
                
                %Get this subject's surface and volume task data 
                %surfdata = [datadir '/' subjects{s} '_mergeFile_L_32k_fsLR.func.gii'];
                surfdata = [datadir '/' subjects{s} '_Mergefile_L.func.gii'];
                voldata = [datadir '/' subjects{s} '_Mergefile.4dfp.img'];
                
                %Track progress
                string{s} = ['    Division ' num2str(divisionnum) ' of ' num2str(divisions) ', subject #' num2str(s) ': ' subjects{s}];
                if s==1; fprintf('%s',string{s}); else fprintf([repmat('\b',1,length(string{s-1})) '%s'],string{s}); end
                
                %Convert the surface to be matlab-readable and the volume to 4dfp
                if divisionnum == 1
                    evalc(['!caret_command64 -file-convert -format-convert ASCII ' surfdata]);
                    delete([surfdata(1:end-9) '_noHEAD.func.gii']);
                    evalc(['!awk ''NF > 25'' ' surfdata ' > ' surfdata(1:end-9) '_noHEAD.func.gii']);
                    
                end
                
                %Read in the surface data and get the betas to be parcellated
                surf_BOLD = load([surfdata(1:end-9) '_noHEAD.func.gii']);
                surf_BOLD(:,1) = [];
                %surf_BOLD = surf_BOLD(:,betastouse+1);
                
                %Read in the volume data and get the betas to be parcellated 
                %voldata4dfp = [voldata(1:end-7) '.4dfp.img'];
                BOLD = read_4dfpimg(voldata);
                vol_BOLD = BOLD(logical(volmask),:)';
                vol_BOLD(find(vol_BOLD==0)) = NaN;
                
                %Set up the avgcrosscorr variable (avg across subjects) for this division
                if s==1
                    avgcrosscorr = zeros(length(indices{divisionnum}), size(vol_BOLD,2));
                end
                
                %Calculate the surface-volume correlation
                crosscorr = FisherTransform(paircorr_mod(single(surf_BOLD(indices{divisionnum},:))',single(vol_BOLD)));
                
                %add that correlation to the avgcrosscorr variable
                avgcrosscorr = avgcrosscorr + crosscorr;
                
                clear crosscorr
                
            end
            
            %Track progress
            string{end+1} = ['    Division ' num2str(divisionnum) ' of ' num2str(divisions) ' complete'];
            fprintf([repmat('\b',1,length(string{s})) '%s'],string{end});
            clear string
            disp(' ')
            
            %Divide the avgcrosscorr by the number of subjects
            avgcrosscorr = avgcrosscorr ./ length(subjects);
            
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
        %delete *mergeFile.nii
        %delete *mergeFile.4dfp*
        delete *_noHEAD.func.gii
        
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
        system(['rm ' outputdir '/CrosscorrDivision_*.mat']);
        
%         avgcorrofcorr = paircorr_mod(single(avgcrosscorr)');
%         clear avgcrosscorr
        cd(outputdir)
        
%         avgcorrofcorrStevedata = avgcorrofcorr;
%         save([outputdir '/avgcorrofcorrStevedata.mat'],'avgcorrofcorrStevedata','-v7.3')
%         clear all
%         surface_taskparcellation_bysubs
        
        
        avgcorrofcorr(isnan(avgcorrofcorr)) = 0;
%        avgcorrofcorr(medialind,:) = [];
        
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
