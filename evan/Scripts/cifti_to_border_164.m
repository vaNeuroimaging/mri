function cifti_to_border_164(ciftiname)
%cifti_to_border(ciftiname)



data_init = ft_read_cifti_mod(ciftiname); data_init = data_init.data;
maskL = gifti(['/data/cn4/evan/Conte_164k/L.atlasroi.164k_fs_LR.func.gii']); maskL = maskL.cdata;
maskR = gifti(['/data/cn4/evan/Conte_164k/R.atlasroi.164k_fs_LR.func.gii']); maskR = maskR.cdata;
data = zeros(length(maskL)+length(maskR) , size(data_init,2));
data(logical([maskL;maskR]),:) = data_init(1:nnz([maskL;maskR]),:);
clear data_init

% if ~exist('bordercolor')
%     bordercolor = [0 0 0];
% end

    


extensionloc = strfind(ciftiname,'.dtseries.nii');

hems = {'L','R'};
hemcaps = {'LEFT', 'RIGHT'};


bufsize=16384;
caretdir = '/data/cn4/evan/Conte_164k/';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

for hemnum = 1:length(hems)
    
    thismask = gifti(['/data/cn4/evan/Conte_164k/' hems{hemnum} '.atlasroi.164k_fs_LR.func.gii']); thismask = thismask.cdata;

outname = [ciftiname(1:(extensionloc-1)) '_' hems{hemnum} '.border'];

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


MNI{hemnum} = gifti(['/data/cn4/evan/Conte_164k/Conte69.' hems{hemnum} '.midthickness.164k_fs_LR.coord.gii']); MNI{hemnum} = MNI{hemnum}.vertices;

dataverts = [1:length(thismask)] + (length(maskL)*(hemnum-1));

for col = 1:size(data,2)
    stufftowritehere = {'   <Border>',...       
        ['      <Name>' num2str(col) '</Name>'],...
        ['      <ClassName>' ciftiname(1:(extensionloc-1)) '</ClassName>'],...
        '      <ColorName>CLASS</ColorName>'};
    for i = 1:length(stufftowritehere)
        dlmwrite(outname,stufftowritehere{i},'-append','delimiter','');
    end
    

    
    verts = find(data(dataverts,col));
    for vert = verts(:)'
        
        %hem = (vert > length(maskL)) + 1;
        
        %vert_withinhem = vert - (length(maskL) * (hem-1));
                
        thirdneighbors = intersect(neighbors(vert,3:end),neighbors(neighbors(vert,2),2:end));
        
        stufftowritehere = {'       <SurfaceProjectedItem>',...
            ['          <Structure>CORTEX_' hemcaps{hemnum} '</Structure>'],...
            ['          <StereotaxicXYZ>' num2str(MNI{hemnum}(vert,:)) '</StereotaxicXYZ>'],...
            '           <ProjectionBarycentric>',...
            '              <TriangleAreas>1 0 0</TriangleAreas>',...
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
end