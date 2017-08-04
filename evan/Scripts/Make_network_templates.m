tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_TMASKLIST.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

%networks = cifti_read('/data/cn4/evan/RestingState/FC_Mapping_120/120_LR_vertexwise_infomap/120_LR_minsize400_recolored.dtseries.nii');
%networks = cifti_read('/data/cn4/evan/RestingState/FC_Mapping_120/120_LR_vertexwise_infomap/120_LR_minsize400_recolored_manualconsensus.dtseries.nii');
networks = cifti_read('/data/cn4/evan/fsaverage_LR32k/Yeo_17_powercolors.dtseries.nii');

templatefile = '/data/hcp-zfs/shared-nil/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';

for s = 1:length(subjects)
    disp(['Subject ' num2str(s)])
    subject = subjects{s};
    cifti_file = ['/data/hcp-zfs/shared-nil/laumannt/120_parcellation/cifti_timeseries_normalwall/' subject '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
    subdata = cifti_read(cifti_file);
    tmask = load(tmasks{s});
    subdata = subdata(:,logical(tmask));
    
    for thresh = 1:size(networks,2)
        vals = unique(networks(:,thresh)); vals(vals<1) = [];
        
        for valnum = 1:length(vals)
            
            correlation = paircorr_mod(subdata',mean(subdata(networks(:,thresh)==vals(valnum),:),1)');
            correlation(isnan(correlation)) = 0;
            correlation = FisherTransform(correlation);
            
            if s==1
                templates{thresh}(:,valnum) = (correlation / length(subjects));
                IDs{thresh}(valnum) = vals(valnum);
            else
                templates{thresh}(:,valnum) = templates{thresh}(:,valnum) + (correlation / length(subjects));
            end
            
        end
    end
end

save('Templates_Yeo.mat','templates','IDs')

for thresh = 1:size(networks,2)
    cifti_write_wHDR(templates{thresh},templatefile,['Template_Yeo'])
end