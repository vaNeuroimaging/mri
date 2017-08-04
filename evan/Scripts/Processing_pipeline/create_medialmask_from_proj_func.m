function create_medialmask_from_proj_func(tmasklist,funcdir,outputdir)
% Create medial wall mask based on surface projection, TOL 09/14
HEMS = {'L';'R'};
smooth = 2.55;
[sessions tmasks] = textread(tmasklist,'%s%s');
for h = 1:2
    for s = 1:length(sessions)
        
        disp(['Creating medialmask, processing subject #' num2str(s)])
        
        timename = [funcdir '/surf_timecourses/' sessions{s} '_' HEMS{h} '_dil10_32k_fs_LR_smooth' num2str(smooth)];
        
        system(['caret_command64 -file-convert -format-convert ASCII ' timename '.func.gii'])
        system(['rm ' timename '_noHEAD.func.gii']);
        system(['awk ''NF > 25'' ' timename '.func.gii > ' timename '_noHEAD.func.gii'])
          
        % Load surface
        surf_time = load([timename '_noHEAD.func.gii']);
        surf_time(:,1) = [];
        noproj(:,s) = sum(surf_time==0,2)==size(surf_time,2);

    end
    sum_noproj = sum(noproj,2)>0;
    save(gifti(single(sum_noproj)),[outputdir '/medial_exclude_' HEMS{h} '.func.gii'])
    atlas_roi = ~logical(sum_noproj);
    save(gifti(single(atlas_roi)),[outputdir '/' HEMS{h} '.atlasroi_noproj.func.gii'])
end