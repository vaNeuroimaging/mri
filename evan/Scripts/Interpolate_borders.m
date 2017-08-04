function Interpolate_borders(borderfile,numbetween)

addpath /data/cn4/evan/Scripts/XML4MATv2/
warning off

lines = textread(borderfile,'%s','delimiter','\t');
newlines = lines;
nodelines = 7;
amounttomove = nodelines * numbetween;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
%neighbors = neighbors+1;

for line = length(lines) : -1 : 1
    
    if strcmp(lines{line},'</Border>')
        first = 1;
    
    elseif strcmp(lines{line},'<SurfaceProjectedItem>')
        if first==1
            first = 0;
        else
            newlines(line+nodelines+amounttomove : end+amounttomove) = newlines(line+nodelines : end);
            newlines(line+nodelines:line+nodelines+amounttomove-1) = repmat(newlines(line:line+nodelines-1),numbetween,1);
            
            thisnodetext = lines{line+3};
            digitlocs = regexp(thisnodetext,'\d');
            spacelocs = regexp(thisnodetext,' ');
            thisnode = thisnodetext(digitlocs(1) : spacelocs(1)-1);
            
            nextnodetext = lines{line+10};
            digitlocs = regexp(nextnodetext,'\d');
            spacelocs = regexp(nextnodetext,' ');
            nextnode = nextnodetext(digitlocs(1) : spacelocs(1)-1);
            
            thisnode_neighs = neighbors(str2num(thisnode),2:7); thisnode_neighs(thisnode_neighs==str2num(nextnode)) = [];
            thirdnode = num2str(thisnode_neighs(1));
            
            
            for i = 1:numbetween
                newlines{line + 3 + 7*i} = ['<TriangleNodes>' thisnode ' ' nextnode ' ' thirdnode '</TriangleNodes>'];
                %newlines{line + 4 + 7*i} = ['<TriangleAreas>' num2str(i/(numbetween+1)) ' ' num2str(1-(i/(numbetween+1))) ' 0.0</TriangleAreas>'];
                newlines{line + 4 + 7*i} = ['<TriangleAreas>1.0 1.0 0.001</TriangleAreas>'];
            end
        end
    end
end
            
filextensionloc = strfind(borderfile,'.border');
outputfilename = [borderfile(1:(filextensionloc-1)) '_interpolated.border'];

dlmwrite(outputfilename,'');   

for line = 1:length(newlines)
    dlmwrite(outputfilename,newlines{line},'-append','delimiter','');
    if ~isempty(strfind(newlines{line},'<TriangleAreas>'))
        dlmwrite(outputfilename,'<SignedDistanceAboveSurface>0</SignedDistanceAboveSurface>','-append','delimiter','');
    end
end


% disp('Reading XML file')
% [y,varname]=xml2struct(borderfile);
% 
% for bordernum = 1:length(y)
%     disp(['Interpolating within border number ' num2str(bordernum)])
%     thisborder = y(bordernum).Border;
%     for borderpoint = [(length(thisborder)-1) : -1 : 1];
%         
%         
%         thispointvertstr = thisborder(borderpoint).SurfaceProjectedItem.ProjectionBarycentric.TriangleNodes;
%         thispointvert = str2num(thispointvertstr); thispointvert = thispointvert(1);
%         
%         nextpointvertstr = thisborder(borderpoint+1).SurfaceProjectedItem.ProjectionBarycentric.TriangleNodes;
%         nextpointvert = str2num(nextpointvertstr); nextpointvert = nextpointvert(1);
%         
%         thisborder((borderpoint+numbetween+1) : (end+numbetween)) = thisborder(borderpoint+1:end);
%                 
%         for between = 1:numbetween
%         
%             thisborder(borderpoint+between) = thisborder(borderpoint);
%             thisborder(borderpoint+between).SurfaceProjectedItem.ProjectionBarycentric.TriangleNodes = num2str([thispointvertstr nextpointvertstr nextpointvertstr]);
%             thisborder(borderpoint+between).SurfaceProjectedItem.ProjectionBarycentric.TriangleAreas = num2str([between/(numbetween+1) 1-(between/(numbetween+1)) 0]);
%         
%         end
%     end
%     
%     thispointvertstr = thisborder(end).SurfaceProjectedItem.ProjectionBarycentric.TriangleNodes;
%     thispointvert = str2num(thispointvertstr); thispointvert = thispointvert(1);
%     
%     nextpointvertstr = thisborder(1).SurfaceProjectedItem.ProjectionBarycentric.TriangleNodes;
%     nextpointvert = str2num(nextpointvertstr); nextpointvert = nextpointvert(1);
%     
%     for between = 1:numbetween
%         
%         thisborder(end+1) = thisborder(end);
%         thisborder(end).SurfaceProjectedItem.ProjectionBarycentric.TriangleNodes = num2str([thispointvertstr nextpointvertstr nextpointvertstr]);
%         thisborder(end).SurfaceProjectedItem.ProjectionBarycentric.TriangleAreas = num2str([between/(numbetween+1) 1-(between/(numbetween+1)) 0]);
%         
%     end
%     y(bordernum).Border = thisborder;
% end
% 
% disp('Writing out interpolated border file')
% 
% newXML = mat2xml(y,varname);
% 
% taglocs = strfind(newXML,'<');
% endtaglocs = strfind(newXML,'</');
% diffline_endtaglocs = strfind(newXML,'></') + 1;
% sameline_endtaglocs = setdiff(endtaglocs,diffline_endtaglocs);
% begintaglocs = union(setdiff(taglocs,sameline_endtaglocs),diffline_endtaglocs);
% 
% filextensionloc = strfind(borderfile,'.border');
% outputfilename = [borderfile(1:(filextensionloc-1)) '_interpolated.border'];
% 
% dlmwrite(outputfilename,'');
% 
% for tagnum = 1:length(begintaglocs)
%     if tagnum<length(begintaglocs)
%         tagtext = newXML(begintaglocs(tagnum) : (begintaglocs(tagnum+1)-1));
%     else
%         tagtext = newXML(begintaglocs(tagnum) : end);
%     end
%     
%     spaceloc = strfind(tagtext,' ');
%     endloc = strfind(tagtext,'>');
%     if numel(spaceloc)>0
%         try
%         tagtext(spaceloc(1) : (endloc(1)-1)) = [];
%         catch
%             1;
%         end
%     end
%     
%     tagtext = regexprep(tagtext,'#46;','.');
%     tagtext = regexprep(tagtext,'#32;',' ');
%     tagtext = regexprep(tagtext,'#95;','_');
%     
%     dlmwrite(outputfilename,tagtext,'-append','delimiter','');
%     
% end


    