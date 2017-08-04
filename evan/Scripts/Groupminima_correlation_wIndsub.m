
groupminimafile = '/data/cn4/evan/RestingState/FC_Mapping_120/MinimaAssignments/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_minima_dist3.func.gii';
groupconnectivityfile = '/data/cn4/evan/RestingState/FC_Mapping_120/MinimaAssignments/avgcorrelmat.mat';
tmaskname = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_DATALIST_TMASK.txt';
medialmask = '/data/cn4/evan/RestingState/FC_Mapping_120/medialmask.func.gii';
groupedgemapfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/';
cohortfile = '/data/cn4/evan/RestingState/FC_Mapping_120/C1toC3_surfaceparcellation.list';

outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/MinimaAssignments/';
outputstatsfile = [outputfolder '/Minimacorrelations.txt'];

delete(outputstatsfile);
fid = fopen(outputstatsfile,'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\n\r\','Subject','GroupNode','Correlation'); %write the output file header
fclose(fid);
dlmwrite(outputstatsfile,' ','-append');

surfacecoordfile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.coord.gii';
surfacetopofile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.32k_fs_LR.topo.gii';

medialmaskdata = gifti(medialmask);
medialmaskind = find(~medialmaskdata.cdata);

groupminimadata = gifti(groupminimafile);
groupminimadata = groupminimadata.cdata;
groupminimadata(medialmaskind) = 0;
groupminimavertices = find(groupminimadata);

load(groupconnectivityfile)
avgcorrelmat(:,medialmaskind) = [];

[tmasksubjects tmaskfiles]=textread(tmaskname,'%s%s');

[ign subnum] = system(['cat ' cohortfile ' | wc -l']);
subnum = str2num(subnum);


disp('Calculating average connectivity patterns of group minima');

for s = 1:subnum
    [ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
    [ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
    [ign funcvols{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1}''']);
    [ign surfdirs{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' cohortfile]);
    
    surfBOLDfilename = dir([groupedgemapfolder '/' subjects{s} '_BOLD*.func.gii']);
    surfBOLDfile{s} = [groupedgemapfolder '/' surfBOLDfilename(1).name(1:end-9)];
    
    tmask{s} = load(tmaskfiles{s});
    
    disp(['    Subject ' subjects{s}])
    evalc(['!rm ' surfBOLDfile{s} '_noHEAD.func.gii']);
    evalc(['!awk ''NF > 25'' ' surfBOLDfile{s} '.func.gii > ' surfBOLDfile{s} '_noHEAD.func.gii']);
    surf_BOLD = load([surfBOLDfile{s} '_noHEAD.func.gii']);
    surf_BOLD(:,1) = [];
    surf_BOLD = single(surf_BOLD(:,logical(tmask{s})));
    evalc(['!rm ' surfBOLDfile{s} '_noHEAD.func.gii']);
    
    subjcorrelmat = FisherTransform(paircorr_mod(surf_BOLD(groupminimavertices,:)',surf_BOLD'));
    subjcorrelmat(:,medialmaskind) = [];
    
    
    %groupvssubcorrelmat = paircorr_mod(subjcorrelmat,avgcorrelmat);
    
    for minimanum = 1:length(groupminimavertices)
        
        correlval = corr2(subjcorrelmat(minimanum,:),avgcorrelmat(minimanum,:));
        
        texttowrite = [subjects{s} '   ' num2str(groupminimavertices(minimanum)) '   ' num2str(correlval)];
        dlmwrite([outputstatsfile],texttowrite,'-append','delimiter','');
    end
    
    
    
end


