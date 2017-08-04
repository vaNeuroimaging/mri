function metric_to_border(metricname,hem)
%metric_to_border(metricname,hem)

metric = gifti(metricname); metric = metric.cdata;

extensionloc = strfind(metricname,'.func.gii');

outname = [metricname(1:(extensionloc-1)) '.border'];


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
    '      <Label Key="0" Red="1" Green="1" Blue="1" Alpha="0"><![CDATA[???]]></Label>',...
    '   </LabelTable>'};

for i = 1:length(header)
    dlmwrite(outname,header{i},'-append','delimiter','');
end

if strcmp(hem,'L') || strcmp(hem,'left') || strcmp(hem,'LEFT')
    hemcaps = 'LEFT';
    hemshort = 'L';
elseif strcmp(hem,'R') || strcmp(hem,'right') || strcmp(hem,'RIGHT')
    hemcaps = 'RIGHT';
    hemshort = 'R';
end

MNI = gifti(['/data/cn4/evan/fsaverage_LR32k/Conte69.' hemshort '.midthickness.32k_fs_LR.surf.gii']); MNI = MNI.vertices;

for col = 1:size(metric,2)
    stufftowritehere = {'   <Border>',...
        ['      <Name>' num2str(col) '</Name>'],...
        '      <ClassName>???</ClassName>',...
        '      <ColorName>CLASS</ColorName>'};
    for i = 1:length(stufftowritehere)
        dlmwrite(outname,stufftowritehere{i},'-append','delimiter','');
    end
    

    
    verts = find(metric(:,col));
    for vert = verts(:)'
        
        thirdneighbors = intersect(neighbors(vert,3:end),neighbors(neighbors(vert,2),2:end));
        
        stufftowritehere = {'       <SurfaceProjectedItem>',...
            ['          <Structure>CORTEX_' hemcaps '</Structure>'],...
            ['          <StereotaxicXYZ>' num2str(MNI(vert,:)) '</StereotaxicXYZ>'],...
            '           <ProjectionBarycentric>',...
            '              <TriangleAreas>1 .001 .001</TriangleAreas>',...
            ['              <TriangleNodes>' num2str([vert-1 neighbors(vert,2)-1 thirdneighbors(1)-1]) '</TriangleNodes>'],...
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