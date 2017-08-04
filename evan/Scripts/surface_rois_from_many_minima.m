minimafolder = '/data/cn4/evan/RestingState/FC_Mapping_120/';
maskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT.func.gii';
outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/Old/';
outputfilename = 'allminima_4mm';
ROIradius = 4;
Hem = 'L';


minimafiles = dir([minimafolder '/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_minima4.func.gii']);
for minimanum = 1:length(minimafiles)
    thisminima = gifti([minimafolder '/' minimafiles(minimanum).name]);
    allminimamatrix(:,minimanum) = thisminima.cdata;
end

minimadata = double(logical(sum(allminimamatrix,2)));


if ~isempty(maskfile)

maskdata = gifti(maskfile);
maskdata = maskdata.cdata;

minimadata = minimadata .* (~maskdata);

end


delete([outputfolder 'minimalocations.txt']);
fid = fopen([outputfolder 'minimalocations.txt'],'at');
fclose(fid);
dlmwrite([outputfolder 'minimalocations.txt'],'','-append');

delete([outputfolder 'minimanames.txt']);
fid = fopen([outputfolder 'minimanames.txt'],'at');
fclose(fid);
dlmwrite([outputfolder 'minimanames.txt'],'','-append');

networkcounter = [];

minimas = find(minimadata)';

for minimaloc = minimas(1:2)

    %dlmwrite([outputfolder 'minimalocations.txt'],num2str(minimaloc),'-append','delimiter','');
    
    orig_minimanum = find(allminimamatrix(minimaloc,:));
    orig_minimanum = orig_minimanum(1);
    
    
    if length(networkcounter) < orig_minimanum
        networkcounter(orig_minimanum) = 0;
    end
    networkcounter(orig_minimanum) = networkcounter(orig_minimanum) +1;
    
    minimaname = [minimafiles(orig_minimanum).name(1:end-9) '_' num2str(networkcounter(orig_minimanum))];
    
    dlmwrite([outputfolder '/minimanames.txt'],minimaname,'-append','delimiter','');

end
dlmwrite([outputfolder '/minimalocations.txt'],num2str(minimas(1:2) -1),'-append','delimiter','');

system(['/data/cn4/evan/workbench/bin_linux64/wb_command -surface-geodesic-rois /data/cn4/evan/fsaverage_LR32k/Conte69.' Hem '.midthickness.32k_fs_LR.surf.gii ' num2str(ROIradius) ' ' outputfolder '/minimalocations.txt ' outputfolder '/' outputfilename '_rois.' Hem '.func.gii -names ' outputfolder 'minimanames.txt -overlap-logic CLOSEST']);
    
    
    
    