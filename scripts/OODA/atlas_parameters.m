function atlas_params = atlas_parameters(atlas)
%
% Name:atlas_parameters.m
% $Revision: 1.3 $
% $Date: 2015/08/04 17:01:15 $
%
homedir = '/data/cn/data1/scripts/OODA/AtlasParameters';
if strcmp(atlas,'Power')
    atlas_params.atlas = atlas;
    atlas_params.networks = {'Unlabelled';'SM';'SM (lat)';'CO';'Auditory';'DMN';'Memory';'Visual';'FP';'Salience';'Sub-cortex';'Ventral Attn';'Dorsal Attn';'Cerebellum'};
    atlas_params.mods = {1:28;29:58;59:63;64:77;78:90;91:148;149:153;154:184;185:209;210:227;228:240;241:249;250:260;261:264};
    atlas_params.atlas_file = 'BigBrain264TimOrder_roilist';
    atlas_params.roi_file = [homedir  '/' atlas_params.atlas_file '_RMAT_fast.roi'];
    load([homedir '/Power_module_colors.mat']); %colors_new
    atlas_params.colors = colors_new;
    
    atlas_params.num_rois = 264;
    atlas_params.dist_thresh = 20;
    
    atlas_params.transitions = [29 59 64 78 91 149 154 185 210 228 241 250 261];
    atlas_params.centers = [14 44 61 70 84 119 151 168 197 218 234 244 255 262];
    atlas_params.sorti = 1:264;

elseif strcmp(atlas,'Parcels')
    atlas_params.atlas = atlas;
    atlas_params.atlas_file = 'Parcels_711-2b_withcoords_roilist';
    atlas_params.roi_file = [homedir '/parcel_center_7112b.roi'];
    atlas_params.dist_thresh = 20; %Tim & Evan typically use "30" for surface distances
    atlas_params.num_rois = 333;
    
    colors_new = [.67 .67 .67;1 0 0;0 0 .6;1 1 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;.8 .8 .6];
    assignments = load([homedir '/ParcelCommunities.txt']);
    IDs = unique(assignments);    
    atlas_params.colors = colors_new(IDs+1,:);
    
    for i = 1:length(IDs)
        mods{i} = find(assignments == IDs(i));
    end
    atlas_params.mods = mods;

    networklabels = {'Unassigned','Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MotorHand','MotorMouth','Auditory','MTL1','MTL2','PERN','RetroSpl'};
    atlas_params.networks = networklabels(IDs+1);
    
    [communities atlas_params.sorti] = sort(assignments);
    atlas_params.transitions = find(communities(1:end-1) ~= communities(2:end));
    transitions_plusends = [1 atlas_params.transitions(:)' length(communities)];
    atlas_params.centers = transitions_plusends(1:end-1) + ((transitions_plusends(2:end) - transitions_plusends(1:end-1))/2);

elseif strcmp(atlas,'ParcelCenters')
    atlas_params.atlas = atlas;
    atlas_params.atlas_file = 'ParcelCenters_roilist';
    atlas_params.roi_file = [homedir '/ParcelCenters_roilist_rmat_slow.roi'];
    atlas_params.dist_thresh = 20; %Tim & Evan typically use "30" for surface distances
    atlas_params.num_rois = 333;
    
    colors_new = [.67 .67 .67;1 0 0;0 0 .6;1 1 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;.8 .8 .6];
    assignments = load([homedir '/ParcelCommunities.txt']);
    IDs = unique(assignments);    
    atlas_params.colors = colors_new(IDs+1,:);
    
    for i = 1:length(IDs)
        mods{i} = find(assignments == IDs(i));
    end
    atlas_params.mods = mods;

    networklabels = {'Unassigned','Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MotorHand','MotorMouth','Auditory','MTL1','MTL2','PERN','RetroSpl'};
    atlas_params.networks = networklabels(IDs+1);
    
    [communities atlas_params.sorti] = sort(assignments);
    atlas_params.transitions = find(communities(1:end-1) ~= communities(2:end));
    transitions_plusends = [1 atlas_params.transitions(:)' length(communities)];
    atlas_params.centers = transitions_plusends(1:end-1) + ((transitions_plusends(2:end) - transitions_plusends(1:end-1))/2);
end

end
