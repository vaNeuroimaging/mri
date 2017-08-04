function cifti_to_border(ciftiname)
%cifti_to_border(ciftiname)



data_init = cifti_read(ciftiname);
maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = ~maskL.cdata;
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = ~maskR.cdata;
data = zeros(length(maskL)+length(maskR) , size(data_init,2));
data(logical([maskL;maskR]),:) = data_init(1:nnz([maskL;maskR]),:);
clear data_init

% if ~exist('bordercolor')
%     bordercolor = [0 0 0];
% end

    


extensionloc = strfind(ciftiname,'.dtseries.nii');

outname = [ciftiname(1:(extensionloc-1)) '.border'];


bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

delete(outname)
fid = fopen(outname,'at'); %open the output file for writing
fclose(fid);

header = {'<?xml version="1.0" encoding="UTF-8"?>',...
    '<BorderFile Version="1">',...
    '   <MetaData>',...
    '      <MD>',...
    '         <Name><![CDATA[UniqueID]]></Name>',...
    '         <Value><![CDATA[{f0f2bfce-760e-4474-b8e4-2a3b1e963eff}]]></Value>',...
    '      </MD>',...
    '   </MetaData>',...
    '   <LabelTable>',...
    ['      <Label Key="1" Red="1" Green="1" Blue="1" Alpha="1"><![CDATA[' ciftiname(1:(extensionloc-1)) ']]></Label>'],...
    '   </LabelTable>'};

for i = 1:length(header)
    dlmwrite(outname,header{i},'-append','delimiter','');
end


hemcaps = {'LEFT', 'RIGHT'};

MNI{1} = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii']); MNI{1} = MNI{1}.vertices;
MNI{2} = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii']); MNI{2} = MNI{2}.vertices;

for col = 1:size(data,2)
    stufftowritehere = {'   <Border>',...       
        ['      <Name>' num2str(col) '</Name>'],...
        ['      <ClassName>' ciftiname(1:(extensionloc-1)) '</ClassName>'],...
        '      <ColorName>CLASS</ColorName>'};
    for i = 1:length(stufftowritehere)
        dlmwrite(outname,stufftowritehere{i},'-append','delimiter','');
    end
    

    
    verts = find(data(:,col));
    for vert = verts(:)'
        
        hem = (vert > length(maskL)) + 1;
        
        vert_withinhem = vert - (length(maskL) * (hem-1));
                
        thirdneighbors = intersect(neighbors(vert_withinhem,3:end),neighbors(neighbors(vert_withinhem,2),2:end));
        
        stufftowritehere = {'       <SurfaceProjectedItem>',...
            ['          <Structure>CORTEX_' hemcaps{hem} '</Structure>'],...
            ['          <StereotaxicXYZ>' num2str(MNI{hem}(vert_withinhem,:)) '</StereotaxicXYZ>'],...
            '           <ProjectionBarycentric>',...
            '              <TriangleAreas>1 0 0</TriangleAreas>',...
            ['              <TriangleNodes>' num2str([vert_withinhem-1 neighbors(vert_withinhem,2)-1 thirdneighbors(1)-1]) '</TriangleNodes>'],...
            '              <SignedDistanceAboveSurface>0</SignedDistanceAboveSurface>',...
            '           </ProjectionBarycentric>',...
            '       </SurfaceProjectedItem>'};
        
        for i = 1:length(stufftowritehere)
            dlmwrite(outname,stufftowritehere{i},'-append','delimiter','');
        end
    end
    dlmwrite(outname,'   </Border>','-append','delimiter','');
end
dlmwrite(outname,'</BorderFile>','-append','delimiter','');