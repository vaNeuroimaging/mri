minimafile = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_minima.func.gii';
maskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT.func.gii';
labelsfile = '/data/cn4/evan/ROIs/FinalLabels/Power_Neuron11_dil.L.32k_fs_LR.label.gii';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/';
outputfilename = 'minima_1mm';
ROIradius = 1;
Hem = 'L';


minimadata = gifti(minimafile);
minimadata = minimadata.cdata;

if ~isempty(maskfile)

maskdata = gifti(maskfile);
maskdata = maskdata.cdata;

minimadata = minimadata .* maskdata;

end

labelsdata = gifti(labelsfile);
labelsdata = labelsdata.cdata;


delete([outputfolder 'minimalocations.txt']);
fid = fopen([outputfolder 'minimalocations.txt'],'at');
fclose(fid);
dlmwrite([outputfolder 'minimalocations.txt'],'','-append');

delete([outputfolder 'minimanames.txt']);
fid = fopen([outputfolder 'minimanames.txt'],'at');
fclose(fid);
dlmwrite([outputfolder 'minimanames.txt'],'','-append');

networkcounter = [];

for minimaloc = find(minimadata)'

    %dlmwrite([outputfolder 'minimalocations.txt'],num2str(minimaloc),'-append','delimiter','');
    
    networknum = labelsdata(minimaloc);
    
    if length(networkcounter) < networknum
        networkcounter(networknum) = 0;
    end
    networkcounter(networknum) = networkcounter(networknum) +1;
    
    [ign networkstring] = system(['awk ''/<Label Key="' num2str(networknum) '"/{print}'' ' labelsfile]);
    stringindex1 = max([strfind(networkstring,['<![CDATA[u' num2str(networknum)]) strfind(networkstring,['<![CDATA[a' num2str(networknum)])]) + length(['<![CDATA[u' num2str(networknum) '_']);
    stringindex2 = strfind(networkstring,']]') - 1;
    networkname = [networkstring(stringindex1:stringindex2) '_' num2str(networkcounter(networknum))];
    
    dlmwrite([outputfolder 'minimanames.txt'],networkname,'-append','delimiter','');

end
dlmwrite([outputfolder 'minimalocations.txt'],num2str(find(minimadata)' -1),'-append','delimiter','');

system(['/data/cn4/evan/workbench/bin_linux64/wb_command -surface-geodesic-rois /data/cn4/evan/fsaverage_LR32k/Conte69.' Hem '.midthickness.32k_fs_LR.surf.gii ' num2str(ROIradius) ' ' outputfolder '/minimalocations.txt ' outputfolder '/' outputfilename '_rois.' Hem '.func.gii -names ' outputfolder 'minimanames.txt -overlap-logic CLOSEST']);
    
    
    
    