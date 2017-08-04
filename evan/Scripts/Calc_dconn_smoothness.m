MSCnums = 1:10;
for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    disp(MSCname)
    leftsurf = [];
    rightsurf = [];
    tmasklist = [];
    
    [subjects, tmasks] = textread(tmasklist,'%s%s');
    thissesssmoothness = zeros(length(subjects),1);
    for s = 1:length(subjects)
        tmask = load(tmasks{s});
        datafile = [];
        data = ft_read_cifti_mod(datafile);
        corr = paircorr_mod(data.data(:,logical(tmask))');
        data.data = corr;
        data.dimord = 'pos_pos';
        ft_write_cifti_mod('temp',data);
        [~, out] = system(['wb_command -cifti-estimate-fwhm temp.dconn.nii -merged-volume -surface CORTEX_LEFT ' leftsurf ' -surface CORTEX_RIGHT ' rightsurf]);
        leftinds = strfind(out,'CORTEX_LEFT FWHM: ');
        rightinds = strfind(out,'CORTEX_RIGHT FWHM: ');
        voxinds = strfind(out,'Voxels FWHM: ');
        
        leftsmooth = zeros(length(leftinds),1);
        rightsmooth = zeros(length(rightinds),1);
        for i = 1:length(leftinds)
            leftsmooth(i) = str2num(out((leftinds(i)+length('CORTEX_LEFT FWHM: ')) : (rightinds(i) -1)));
            rightsmooth(i) = str2num(out((rightinds(i)+length('CORTEX_RIGHT FWHM: ')) : (voxinds(i) -1)));
        end
        
        thissesssmoothness(s) = mean(mean([leftsmooth rightsmooth],2),1);
    end
    indsess_smoothness(MSCnum) = mean(thissesssmoothness);
    
    
    dconnfile = [];
    [~, out] = system(['wb_command -cifti-estimate-fwhm ' dconnfile ' -merged-volume -surface CORTEX_LEFT ' leftsurf ' -surface CORTEX_RIGHT ' rightsurf]);
    leftinds = strfind(out,'CORTEX_LEFT FWHM: ');
    rightinds = strfind(out,'CORTEX_RIGHT FWHM: ');
    voxinds = strfind(out,'Voxels FWHM: ');
    
    leftsmooth = zeros(length(leftinds),1);
    rightsmooth = zeros(length(rightinds),1);
    for i = 1:length(leftinds)
        leftsmooth(i) = str2num(out((leftinds(i)+length('CORTEX_LEFT FWHM: ')) : (rightinds(i) -1)));
        rightsmooth(i) = str2num(out((rightinds(i)+length('CORTEX_RIGHT FWHM: ')) : (voxinds(i) -1)));
    end
    
    sessavg_smoothness(MSCnum) = mean(mean([leftsmooth rightsmooth],2),1);
end

disp('Group')
dconnfile = [];
leftsurf = [];
rightsurf = [];
[~, out] = system(['wb_command -cifti-estimate-fwhm ' dconnfile ' -merged-volume -surface CORTEX_LEFT ' leftsurf ' -surface CORTEX_RIGHT ' rightsurf]);
leftinds = strfind(out,'CORTEX_LEFT FWHM: ');
rightinds = strfind(out,'CORTEX_RIGHT FWHM: ');
voxinds = strfind(out,'Voxels FWHM: ');

leftsmooth = zeros(length(leftinds),1);
rightsmooth = zeros(length(rightinds),1);
for i = 1:length(leftinds)
    leftsmooth(i) = str2num(out((leftinds(i)+length('CORTEX_LEFT FWHM: ')) : (rightinds(i) -1)));
    rightsmooth(i) = str2num(out((rightinds(i)+length('CORTEX_RIGHT FWHM: ')) : (voxinds(i) -1)));
end

groupsmoothness = mean(mean([leftsmooth rightsmooth],2),1);

figure;
plotSpread([indsess_smoothness(:);sessavg_smoothness(:);groupsmoothness],'categoryIdx',[zeros(length(MSCnums),1);ones(length(MSCnums),1);2],'categoryColors',{'b','g','r'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)