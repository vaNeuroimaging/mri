function focifilemaker_workbench(nodes,RGB,fociname,classname,focifilename,hem)
%function focifilemaker_workbench(nodes,RGB,fociname,classname,focifilename,hem)
%
%nodes - a list of cortical nodes that should become foci centers
%RGB - a #nodes X 3 RGB color mapping
%fociname - a #nodes long cell array containing strings of foci names. Foci names may not be repeated.
%classname - a #nodes long cell array containing strings of foci classes. Foci classes may be repeated.
%focifilename - the name the foci file will be written to
%hem - 'L' or 'R'

if strcmp(hem,'L')
    fullhemname = 'LEFT';
elseif strcmp(hem,'R')
    fullhemname = 'RIGHT';
end

surface = gifti(['/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
foci = surface.vertices(nodes,:);

fid=fopen(focifilename, 'wt');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<FociFile Version="1">\n');
fprintf(fid,'   <MetaData>\n');
fprintf(fid,'      <MD>\n');
fprintf(fid,'         <Name><![CDATA[Caret-Version]]></Name>\n');
fprintf(fid,'         <Value><![CDATA[5.64]]></Value>\n');
fprintf(fid,'      </MD>\n');
fprintf(fid,'      <MD>\n');
fprintf(fid,'         <Name><![CDATA[Date]]></Name>\n');
fprintf(fid,'         <Value><![CDATA[2012-05-20T20:23:45]]></Value>\n');
fprintf(fid,'      </MD>\n');
fprintf(fid,'      <MD>\n');
fprintf(fid,'         <Name><![CDATA[UniqueID]]></Name>\n');
fprintf(fid,'         <Value><![CDATA[{1f72c4fa-414b-4e75-b310-4d298ae816bd}]]></Value>\n');
fprintf(fid,'      </MD>\n');
fprintf(fid,'      <MD>\n');
fprintf(fid,'         <Name><![CDATA[UserName]]></Name>\n');
fprintf(fid,'         <Value><![CDATA[vanessen]]></Value>\n');
fprintf(fid,'      </MD>\n');
fprintf(fid,'   </MetaData>\n');
fprintf(fid,'   <LabelTable>\n');
fprintf(fid,'      <Label Key="0" Red="0.666667" Green="0.666667" Blue="0.666667" Alpha="0"><![CDATA[???]]></Label>\n');

for f = 1:size(foci,1)
    
    fprintf(fid,'      <Label Key="%i" Red="%4.2f" Green="%4.2f" Blue="%4.2f"><![CDATA[%s]]></Label>\n',f,RGB(f,1),RGB(f,2),RGB(f,3),fociname{f});
end

fprintf(fid,'   </LabelTable>\n');

for f = 1:size(foci,1)
    
      fprintf(fid,'   <Focus Index="%d">\n',f);
      fprintf(fid,'      <Area><![CDATA[]]></Area>\n');
      fprintf(fid,'      <ClassName><![CDATA[%s]]></ClassName>\n',classname{f});
      fprintf(fid,'      <Comment><![CDATA[]]></Comment>\n');
      fprintf(fid,'      <Extent>0</Extent>\n');
      fprintf(fid,'      <Geography><![CDATA[]]></Geography>\n');
      fprintf(fid,'      <Name><![CDATA[%s]]></Name>\n',fociname{f});
      fprintf(fid,'      <RegionOfInterest><![CDATA[]]></RegionOfInterest>\n');
      fprintf(fid,'      <SearchXYZ>%4.2f %4.2f %4.2f</SearchXYZ>\n',foci(f,1),foci(f,2),foci(f,3));
      fprintf(fid,'      <Statistic><![CDATA[]]></Statistic>\n');
      fprintf(fid,'      <SumsIDNumber><![CDATA[-1]]></SumsIDNumber>\n');
      fprintf(fid,'      <SumsRepeatNumber><![CDATA[-1]]></SumsRepeatNumber>\n');
      fprintf(fid,'      <SumsParentFocusBaseID><![CDATA[-1]]></SumsParentFocusBaseID>\n');
      fprintf(fid,'      <SumsVersionNumber><![CDATA[-1]]></SumsVersionNumber>\n');
      fprintf(fid,'      <SumsMSLID><![CDATA[-1]]></SumsMSLID>\n');
      fprintf(fid,'      <AttributeID><![CDATA[-1]]></AttributeID>\n');
      fprintf(fid,'      <StudyMetaDataLinkSet>\n');
      fprintf(fid,'      </StudyMetaDataLinkSet>\n');
      fprintf(fid,'      <SurfaceProjectedItem>\n');
      fprintf(fid,['         <Structure>CORTEX_' fullhemname '</Structure>\n']);
      fprintf(fid,'         <StereotaxicXYZ>%4.2f %4.2f %4.2f</StereotaxicXYZ>\n',foci(f,1),foci(f,2),foci(f,3));
      fprintf(fid,'         <VolumeXYZ>%4.2f %4.2f %4.2f</VolumeXYZ>\n',foci(f,1),foci(f,2),foci(f,3));
      fprintf(fid,'         <ProjectionBarycentric>\n');
      fprintf(fid,'            <TriangleAreas>1 0 0</TriangleAreas>\n');
      fprintf(fid,'            <TriangleNodes>%i %i %i</TriangleNodes>\n',nodes(f)-1,nodes(f)-1,nodes(f)-1);
      fprintf(fid,'            <SignedDistanceAboveSurface>0</SignedDistanceAboveSurface>\n');
      fprintf(fid,'         </ProjectionBarycentric>\n');
      fprintf(fid,'      </SurfaceProjectedItem>\n');
      fprintf(fid,'   </Focus>\n');
end
fprintf(fid,'</FociFile>\n');
