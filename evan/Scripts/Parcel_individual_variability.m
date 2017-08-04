tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_TMASKLIST.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

make_group_correlpatterns = 1;

ciftitemplatefile = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
parcelsname = '120_subsurf_LR_nosmooth_watershedmerge_0.4_tweaked_gooddata.dtseries.nii';
%%
parcels = cifti_read(parcelsname);
parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0) = [];

group_correlpatterns = zeros(size(parcels,1),length(parcelIDs));

if make_group_correlpatterns

for s = 1:length(subjects);
    
    tmask = load(tmasks{s});
    timecourse = cifti_read(['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii']);
    %timecourse = timecourse.cdata;
    timecourse = timecourse(:,logical(tmask));
    timecourse(isnan(timecourse)) = 0;
    
    for parcelnum = 1:length(parcelIDs)
        parcelindices = find(parcels==parcelIDs(parcelnum));
        
        parceltimecourse = mean(timecourse(parcelindices,:),1);
        
        parcel_correlpattern = FisherTransform(paircorr_mod(timecourse',parceltimecourse'));
        parcel_correlpattern(isnan(parcel_correlpattern)) = 0;
        group_correlpatterns(:,parcelnum) = group_correlpatterns(:,parcelnum) + (parcel_correlpattern / length(subjects));
    end

end
    cifti_write_wHDR(group_correlpatterns,ciftitemplatefile,'All_parcel_correlpatterns_gooddata.dtseries.nii');
else
    group_correlpatterns = cifti_read('All_parcel_correlpatterns_gooddata.dtseries.nii');
end
%%
correloutput = zeros(length(parcels),length(subjects));
%avgcorreloutput = zeros(size(parcels));

i = 0;
for s = 1:length(subjects);
    
    tmask = load(tmasks{s});
    timecourse = cifti_read(['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii']);
    %timecourse = timecourse.cdata;
    timecourse = timecourse(:,logical(tmask));
    timecourse(isnan(timecourse)) = 0;
    
    for parcelnum = 1:length(parcelIDs)
        i = i+1;
        string{i} = ['Subject ' num2str(s) ' of ' num2str(length(subjects)) ', parcel ' num2str(parcelnum) ' of ' num2str(length(parcelIDs))];
        if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
        
        parcelindices = find(parcels==parcelIDs(parcelnum));
        
        parceltimecourse = mean(timecourse(parcelindices,:),1);
        
        parcel_correlpattern = FisherTransform(paircorr_mod(timecourse',parceltimecourse'));
        parcel_correlpattern(isnan(parcel_correlpattern)) = 0;
        
        similarity_to_group = paircorr_mod(parcel_correlpattern,group_correlpatterns(:,parcelnum));
        
        correloutput(parcelindices,s) = similarity_to_group;
        
    end
end



cifti_write_wHDR(correloutput,ciftitemplatefile,'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/Subject_similarity_togroup_byparcels')
cifti_write_wHDR(mean(correloutput,2),ciftitemplatefile,'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/Avg_Subject_similarity_togroup_byparcels')
cifti_write_wHDR(std(FisherTransform(correloutput),[],2),ciftitemplatefile,'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/Std_Subject_similarity_togroup_byparcels')


%%

data = cifti_read('Subject_similarity_togroup_byparcels.dtseries.nii');
parcels = cifti_read(parcelsname);
parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_TMASKLIST.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');
for s = 1:length(subjects)
tmask = load(tmasks{s});
timepoints(s) = nnz(tmask);
for parcelnum = 1:length(parcelIDs)
    Zs(parcelnum,s) = mean(data(parcels==parcelIDs(parcelnum),s));
end
end

figure
hold on
plot(timepoints,mean(Zs,1),'k.','MarkerSize',20)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)

%moredatasubs = find(timepoints>=300);
%cifti_write_wHDR(mean(data(:,moredatasubs),2),ciftitemplatefile,'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/Avg_GoodSubject_similarity_togroup_byparcels')

% f = .5;
% randadd = ((rand(size(timepoints))*2) -1) * .00001;
% datain = [(timepoints + randadd)', mean(Zs,1)'];
% evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
% plot(xy(:,1),xy(:,2),'k');
        
    