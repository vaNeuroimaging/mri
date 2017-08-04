minimafile1 = '/data/cn4/evan/RestingState/FC_Mapping_120/smooth255preedge_surfsmooth_ztrans_bandpass_surfreg_32K_AllC_avg_edge_avg_smooth_L_noalone_minima.func.gii';
minimafile2 = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_minima.func.gii';

surfacecoordfile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.coord.gii';
surfacetopofile = '/data/cn4/evan/fsaverage_LR32k/Conte69.L.32k_fs_LR.topo.gii';
%distancesfile = '/data/cn4/evan/fsaverage_LR32k/VertexDistances.L.func.gii';

minimadata1 = gifti(minimafile1);
minimadata1 = minimadata1.cdata;
minimavertices1 = find(minimadata1);

minimadata2 = gifti(minimafile2);
minimadata2 = minimadata2.cdata;
minimavertices2 = find(minimadata2);

outputmetric1 = zeros(size(minimadata1));
outputmetric2 = outputmetric1;

randomorder = randperm(max(length(minimavertices1),length(minimavertices2)));
nomatchnumber = max(length(minimavertices1),length(minimavertices2)) + 50;

for minimavertexnum = 1:length(minimavertices1)
    minimavertex = minimavertices1(minimavertexnum);
    
    displaystring = ['Node ' num2str(minimavertex)];
    if minimavertexnum==1; fprintf('%s',displaystring); else; fprintf([repmat('\b',1,length(['Node ' num2str(minimavertices1(minimavertexnum-1))])) '%s'],displaystring); end
    
    
    outputmetric1(minimavertex) = randomorder(minimavertexnum);

    system(['caret_command64 -surface-geodesic ' surfacecoordfile ' ' surfacetopofile ' temp.func.gii false -node ' num2str(minimavertex-1)]);
    evalc(['!caret_command64 -file-convert -format-convert XML temp.func.gii']);
    distances = gifti('temp.func.gii');
    distances = distances.cdata;
    distancestominima = distances(minimavertices2);
    [mindist, minindex] = min(distancestominima);
    closestminimum = minimavertices2(minindex);
    
    
    system(['caret_command64 -surface-geodesic ' surfacecoordfile ' ' surfacetopofile ' temp.func.gii false -node ' num2str(closestminimum-1)]);
    evalc(['!caret_command64 -file-convert -format-convert XML temp.func.gii']);
    distances = gifti('temp.func.gii');
    distances = distances.cdata;
    distancestominima = distances(minimavertices1);
    [mindist, minindex] = min(distancestominima);
    reciporicalclosestminimum = minimavertices1(minindex);
    
    if reciporicalclosestminimum == minimavertex
        outputmetric2(closestminimum) = randomorder(minimavertexnum);
    else
        outputmetric1(minimavertex) = nomatchnumber;
    end
    
end

for minimavertexnum = 1:length(minimavertices2)
    minimavertex = minimavertices2(minimavertexnum);
    if outputmetric2(minimavertex) == 0;
        outputmetric2(minimavertex) = nomatchnumber;
    end
end

system('rm temp.func.gii');

save(gifti(outputmetric1),[minimafile1(1:end-9) '_matched.func.gii'],'ExternalFileBinary');
save(gifti(outputmetric2),[minimafile2(1:end-9) '_matched.func.gii'],'ExternalFileBinary');
    


