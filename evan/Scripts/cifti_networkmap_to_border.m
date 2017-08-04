function cifti_networkmap_to_border(ciftifile)
%cifti_to_border_v2(ciftiname,filled,allcolors_inonemap,bordercolors)


networklabels = [{'Default','Visual','FrontoPar','Unknown','DorsalAttn','Unknown','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'},repmat({'Unknown'},1,100)];
powercolors = [1 0 0;0 0 .6;.9 .9 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;.85 .85 .85;.5 .5 .3;repmat([.8 .35 .5],100,1)];

data_template = ft_read_cifti_mod(ciftifile);
data = data_template.data;
data(isnan(data)) = 0;

num_wholesurfverts = nnz(data_template.brainstructure <= 2);
data_wholesurf = zeros(num_wholesurfverts,size(data,2));
wholesurf_indswithdata = data_template.brainstructure(1:num_wholesurfverts) > 0;
data_wholesurf(wholesurf_indswithdata,:) = data(1:nnz(wholesurf_indswithdata),:);

data_bothhems = data_wholesurf;

temp = data_bothhems;
IDs = unique(temp); IDs(IDs==0) = [];
data_bothhems = zeros(size(temp,1),length(IDs));
for IDnum = 1:length(IDs)
    ID = IDs(IDnum);
    data_bothhems(:,IDnum) = temp==ID;
    
    if mod(ID,1)==0
        colnames_bothhem{IDnum} = networklabels{ID};
        bordercolors(IDnum,:) = powercolors(ID,:);
    else
        colnames_bothhem{IDnum} = 'Unknown';
        decimal = mod(ID,1);
        bordercolors(IDnum,:) = (powercolors(ceil(ID),:) * decimal) + (powercolors(floor(ID),:) * (1-decimal));
    end
end
  
clear temp




clear neighbors

bufsize=16384;
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread('/home/data/scripts/Resources/node_neighbors.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


hemcaps = {'LEFT', 'RIGHT'};

MNI{1} = gifti(['/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii']); MNI{1} = MNI{1}.vertices;
MNI{2} = gifti(['/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.midthickness.32k_fs_LR.surf.gii']); MNI{2} = MNI{2}.vertices;


for hem = 1:2
    
    data = data_bothhems([1:size(MNI{hem},1)] + (size(MNI{1},1) * (hem-1)),:);
    
    cols_withdata = logical(sum(data,1));
    data = data(:,cols_withdata);
    colnames = colnames_bothhem(cols_withdata);
    thishem_bordercolors = bordercolors(cols_withdata,:);
    
    
    data_temp = zeros(size(data));
    [verts,cols] = find(data);
    for vertnum = 1:size(verts)
        vertneighbors = neighbors(verts(vertnum),:); vertneighbors(isnan(vertneighbors)) = [];
        if any(data(vertneighbors,cols(vertnum))==0)
            data_temp(verts(vertnum),cols(vertnum)) = 1;
        end
    end
    data = data_temp;
    
    
    
    extensionloc = strfind(ciftifile,'.dtseries.nii');
    if isempty(extensionloc)
        extensionloc = strfind(ciftifile,'.dscalar.nii');
    end
    
    outname = [ciftifile(1:(extensionloc-1)) '_' hemcaps{hem}(1) '.border'];
    
    
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
        '   <BorderClassColorTable>',...
        '      <LabelTable>',...
        ['         <Label Key="1" Red="0" Green="0" Blue="0" Alpha="1"><![CDATA[' ciftifile(1:(extensionloc-1)) ']]></Label>'],...
        '      </LabelTable>',...
        '   </BorderClassColorTable>'};
    for i = 1:length(header)
        dlmwrite(outname,header{i},'-append','delimiter','');
    end
    
    
    header = {'   <BorderNameColorTable>',...
        '      <LabelTable>'};
    for col = 1:size(data,2)
        header{end+1} = ['         <Label Key="' num2str(col) '" Red="' num2str(thishem_bordercolors(col,1)) '" Green="' num2str(thishem_bordercolors(col,2)) '" Blue="' num2str(thishem_bordercolors(col,3)) '" Alpha="1"><![CDATA[' colnames{col} ']]></Label>'];
    end
    header(end+1 : end+2) = {'      </LabelTable>',...
        '   </BorderNameColorTable>'};
    
    for i = 1:length(header)
        dlmwrite(outname,header{i},'-append','delimiter','');
    end
    
    
    for col = 1:size(data,2)
        stufftowritehere = {'   <Border>',...
            ['      <Name>' colnames{col} '</Name>'],...
            ['      <ClassName>' ciftifile(1:(extensionloc-1)) '</ClassName>'],...
            '      <ColorName>CLASS</ColorName>'};
        for i = 1:length(stufftowritehere)
            dlmwrite(outname,stufftowritehere{i},'-append','delimiter','');
        end
        
        
        
        verts = find(data(:,col));
        for vert = verts(:)'
            
            thirdneighbors = intersect(neighbors(vert,3:end),neighbors(neighbors(vert,2),2:end));
            
            stufftowritehere = {'       <SurfaceProjectedItem>',...
                ['          <Structure>CORTEX_' hemcaps{hem} '</Structure>'],...
                ['          <StereotaxicXYZ>' num2str(MNI{hem}(vert,:)) '</StereotaxicXYZ>'],...
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