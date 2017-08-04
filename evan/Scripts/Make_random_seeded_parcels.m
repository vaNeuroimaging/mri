parcelname = '/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh2_watershedmerge_0.45.func.gii';
hem = 'R';

medialmaskdata = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']);
medialmaskdata = ~medialmaskdata.cdata;
corticalinds = find(medialmaskdata);

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


parcels = gifti(parcelname); parcels = parcels.cdata;

numwatersheds = nnz(unique(parcels));

label = zeros(32492,1);
%     randorder = randperm(length(corticalinds));
%     randverts = corticalinds(randorder(1:numwatersheds));
%label(randverts) = randverts;
    
for i = 1:numwatersheds
    done = 0;
    while done == 0
        randchoice = corticalinds(randi(length(corticalinds)));
        choiceneighs = neighbors(randchoice,:); choiceneighs(isnan(choiceneighs)) = [];
        if ~any(label(choiceneighs)>0)
            label(randchoice) = randchoice;
            done = 1;
        end
    end
end


    
    stillexpanding = 1;
    while stillexpanding==1
        stillexpanding = 0;
        borderverts = find((label==0) .* medialmaskdata);
        borderverts = borderverts(randperm(length(borderverts)));
        for vert = borderverts'
            vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs)) = [];
            neighlabels = label(vertneighs); neighlabels(neighlabels==0) = [];
            if nnz(unique(neighlabels)) == 1
                stillexpanding = 1;
                label(vert) = neighlabels(1);
            end
        end
    end
    
    save(gifti(single(label)),['/data/cn4/evan/RestingState/FC_Mapping_120/Random_matchedto_120_R_wateredgethresh_watershedmerge_0.45.func.gii'])