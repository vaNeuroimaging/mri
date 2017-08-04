
%dir = '/data/cn4/laumannt/fcMapping_redux';

[subjects tmasks] = textread(['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/60sub_C1_NEW_TMASKLIST.txt'],'%s%s');

surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR';
HEMS = {'L';'R'};

for h = 1:2
    for s = 1:length(subjects)
        
        label = gifti([surfdir '/' subjects{s} '/7112b_fs_LR/fsaverage_LR32k/' subjects{s} '.' HEMS{h} '.aparc.a2009s.32k_fs_LR.label.gii']);
        label = label.cdata;
        labels{h}(:,s) = label;
    end
end

%%
%outputdir = '/data/cn4/laumannt/fcMapping_redux';
for h = 1:2
mode_label = mode(single(labels{h}),2);
save(gifti(single(mode_label)),['C1_mode.' HEMS{h} '.aparc.a2009s.32k_fs_LR.func.gii'])
end