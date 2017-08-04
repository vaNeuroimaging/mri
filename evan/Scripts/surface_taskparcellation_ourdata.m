%function  surface_parcellation(cohortfile,outputdir,tmasktype,smooth,hems,edges)
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

hem = 'R';
hems = 'RIGHT';
outputdir = '/data/cn4/evan/Task_parcellation/Manytasks/';
voldata = '/data/cn4/evan/Task_parcellation/Manytasks/Manytasks.4dfp.img';
surfdata = ['/data/cn4/evan/Task_parcellation/Manytasks/Manytasks_' hem '.func.gii'];
smooth = 'yes';
edges = 'yes';
medialmask = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
medialmaskdata = gifti(medialmask);
medialind = find(medialmaskdata.cdata);

mask4dfp='/data/cn4/laumannt/Standard/glm_atlas_mask_333.4dfp.img';
volmask = read_4dfpimg(mask4dfp);
    

divisions = 12;

smoothnum = 2.55;
%mask4dfp='/data/cn4/laumannt/Standard/glm_atlas_mask_222.4dfp.img';
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};




switch hems
    case 'LEFT'
        h = 1;
    case 'RIGHT'
        h = 2;
    case 'BOTH'
        h = [1:2];
end

% Make output folder
system(['mkdir ' outputdir])

for hem = h
    

    maskout = zeros(1,65549);
    
        cd(outputdir)


        surf_BOLD = gifti(surfdata);
        surf_BOLD = surf_BOLD.cdata;
        
        %mask = read_4dfpimg(mask4dfp);
        BOLD = read_4dfpimg(voldata);
        
        vol_BOLD = BOLD(logical(volmask),:)';        
        
        % Calculate correlation between surface and volume data
        nodesperdivision = ceil(size(surf_BOLD,1) / divisions);
        
        disp('Calculating correlation between surface nodes and volume')
        for divisionnum = 1:divisions
            
            if divisionnum==divisions
                indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : size(surf_BOLD,1);
            else
                indices{divisionnum} = nodesperdivision*(divisionnum-1)+1 : nodesperdivision*divisionnum;
            end
            
            crosscorr = FisherTransform(paircorr_mod(single(surf_BOLD(indices{divisionnum},:))',single(vol_BOLD)));
            
            % Remove voxels from crosscorr that have no BOLD data but are in glm mask
            nanmat = isnan(crosscorr);
            nancolsum = sum(nanmat,2)==size(crosscorr,2);
            nanmat(logical(nancolsum),:) = [];
            nanrowsum = sum(nanmat,1)>0;
            %nanrowsum = repmat(nanrowsum,1,numwindows);
            crosscorr(:,logical(nanrowsum)) = [];
            
            save([outputdir '/CrosscorrDivision' num2str(divisionnum) '.mat'],'crosscorr','-v7.3')
            
            clear crosscorr
            
        end
        
        
        %evalc(['!rm ' BOLDname '_noHEAD.func.gii']);
        %evalc(['!rm ' funcvol '.*']);
        disp('done.')
        
        
        
        %         avgcrosscorr = crosscorr;
        %         clear crosscorr
        
        
        
        
        % Create mask based on nodes that don't get projected to
        %medialind = isnan(medialmaskcheckmatrix);
        %medialmask = ones(length(medialmaskcheckmatrix),1);
        %medialmask(medialind) = 0;
        %save(gifti(single(medialmask)),'medialmask.func.gii')
        
        %for piecenum = 1:(numwindows*2)
        
        % Calculate correlation of correlation maps
        fprintf('Calculating correlation of correlation maps...')
        disp(' ')
        avgcorrofcorr = zeros(size(surf_BOLD,1),'single');
        
        
        
        divisionnum = 0;
        
        for i = 1:divisions
            
            load([outputdir '/CrosscorrDivision' num2str(i) '.mat']);
            crosscorri = crosscorr;
            
            for j = 1:divisions
                
                if i >= j
                    
                    load([outputdir '/CrosscorrDivision' num2str(j) '.mat']);
                    crosscorrj = crosscorr;
                    
                    divisionnum = divisionnum+1;
                    
                    string{divisionnum} = ['Calculating correlations for division number ' num2str(divisionnum) ' of ' num2str((divisions^2 + divisions)/2)];
                    if divisionnum==1; fprintf('%s',string{divisionnum}); else fprintf([repmat('\b',1,length(string{divisionnum-1})) '%s'],string{divisionnum}); end
                    
                    avgcorrofcorr(indices{i},indices{j}) = paircorr_mod(crosscorri',crosscorrj');
                    
                end
                
            end
        end
        disp(' ')
        
        clear crosscorr crosscorri crosscorrj
        
        avgcorrofcorr = tril(avgcorrofcorr,-1) + tril(avgcorrofcorr)';
        
        system(['rm ' outputdir '/CrosscorrDivision*.mat']);
        
        
%         
%         
%         %crosscorr = FisherTransform(paircorr_mod(single(surf_BOLD)',single(vol_BOLD)));
%         
%         %crosscorr_surf = FisherTransform(paircorr_mod(single(surf_BOLD)'));
%         
% %         % Keep running correlation map average
% %         if numel(subjects)>1
% %             avgcrosscorr = avgcrosscorr + crosscorr;
% %             %avgcrosscorr_surf = avgcrosscorr_surf + crosscorr_surf;
% %         end
%         
%         % Keep track of voxels without data by subject
%         maskout = maskout + isnan(crosscorr(1,:));
%         
%         system(['rm ' BOLDname '_noHEAD.func.gii'])
%         system(['rm ' funcvol '.*'])
%         disp('done.')
%     %end
%     
% %     if numel(subjects)>1
% %         clear crosscorr crosscorr_surf
% %         avgcrosscorr = avgcrosscorr./(length(subjects));
% %         %avgcrosscorr_surf = avgcrosscorr_surf./(length(subjects));
% %         %save([filestem 'avgcrosscorr_surf_' HEMS{hem} '.mat'],'avgcrosscorr_surf','-v7.3')
% %     else
%         avgcrosscorr = crosscorr;
%         clear crosscorr
% %    end
%     
%     % Remove voxels from crosscorr that have no BOLD data but are in glm mask
%     nanmat = isnan(avgcrosscorr);
%     nancolsum = sum(nanmat,2)==length(avgcrosscorr(1,:));
%     nanmat(logical(nancolsum),:) = [];
%     nanrowsum = sum(nanmat,1)>0;
%     avgcrosscorr(:,logical(nanrowsum)) = [];
%     
%     % Create mask based on nodes that don't get projected to
%     medialind = isnan(avgcrosscorr(:,1));
%     medialmask = ones(size(avgcrosscorr,1),1);
%     medialmask(medialind) = 0;
%     save(gifti(single(medialmask)),'medialmask.func.gii')
%     
%     % Calculate correlation of correlation maps
%     fprintf('Calculating correlation of correlation maps...')
%     avgcorrofcorr = paircorr_mod(single(avgcrosscorr)');
    disp('done.')
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
    filestem = ['avgcorrofcorr_smooth' num2str(smoothnum) '_allgrad_' HEMS{hem}];
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
    grad_metrics_piece1 = load([outputdir '/' filestem '_piece1_noHEAD.func.gii']);
    grad_metrics_piece1(:,1) = [];
    grad_metrics_piece2 = load([outputdir '/' filestem '_piece2_noHEAD.func.gii']);
    grad_metrics_piece2(:,1) = [];
    grad_metrics = [grad_metrics_piece1 grad_metrics_piece2];
    
    % Calculate and save average gradient map
    avg_metric = mean(grad_metrics,2);
    save(gifti(single(avg_metric)),[filestem '_avg.func.gii']);
    
    % Calculate edge maps
    switch edges
        case 'no'
        case 'yes'
            surface_edges_all_test_faster(grad_metrics,specfile,outputdir,[filestem '_edge_avg'],1);
    end
       
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