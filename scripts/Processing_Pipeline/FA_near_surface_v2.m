function FA_near_surface_v2(subject,distances)


Lpial = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.pial.32k_fs_LR.surf.gii']);
Rpial = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.pial.32k_fs_LR.surf.gii']);
Lwhite = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.white.32k_fs_LR.surf.gii']);
Rwhite = gifti(['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.white.32k_fs_LR.surf.gii']);

bufsize=16384;
surfneighborfile = '/home/data/scripts/Resources/node_neighbors.txt';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[surfneighbors(:,1) surfneighbors(:,2) surfneighbors(:,3) surfneighbors(:,4)...
    surfneighbors(:,5) surfneighbors(:,6) surfneighbors(:,7)] = ...
    textread(surfneighborfile,'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
surfneighbors = surfneighbors+1;

for d = 1:length(distances)
    
    distance = distances(d);
    
    Linteriorsurf = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.interior' num2str(distance) '.32k_fs_LR.surf.gii'];
    Rinteriorsurf = ['/home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.interior' num2str(distance) '.32k_fs_LR.surf.gii'];
    
    %if ~exist(Linteriorsurf,'file') || ~exist(Rinteriorsurf,'file')
        
        vecsL = Lwhite.vertices - Lpial.vertices;
        mags = sqrt(sum((vecsL.^2),2));
        normvecsL = vecsL .* distance ./ repmat(mags,1,3);
        normvecsL(isnan(normvecsL)) = 0;
        interiorL = Lwhite.vertices + normvecsL;
        
        
        
%         white_to_interior_distances = pdist2(Lwhite.vertices,interiorL);
%         toownvert_distances = diag(white_to_interior_distances,0);
%         distance_diffs = repmat(toownvert_distances,1,size(white_to_interior_distances,2)) - white_to_interior_distances;
%         distance_diffs(mags==0,:) = 0;
%         distance_diffs(:,mags==0) = 0;
%         for i = 1:size(distance_diffs,1)
%             theseneighs = surfneighbors(i,:); theseneighs(isnan(theseneighs)) = [];
%             distance_diffs(i,theseneighs) = 0;
%         end
%         smaller_than_toownvert_inds = find(distance_diffs>.01);
%         [~, sortindex] = sort(distance_diffs(smaller_than_toownvert_inds),'ascend');
%         
%         clear distance_diffs white_to_interior_distances
%         
%         for i = 1:length(sortindex)
%             if mod(i,1000)==0
%                 disp([num2str(i) ' of ' num2str(length(sortindex))])
%             end
%             thisind = smaller_than_toownvert_inds(sortindex(i));
%             [vertind, otherind] = ind2sub([length(mags) length(mags)],thisind);
%             whitecoord = Lwhite.vertices(vertind,:);
%             newvertcoord = interiorL(vertind,:);
%             newothercoord = interiorL(otherind,:);
%             
%             for j = [.001 : .001 : 20]
%                 if pdist2(whitecoord,newothercoord) >= pdist2(whitecoord,newvertcoord)
%                     break
%                 else
%                     thisvertvec = vecsL(vertind,:) .* j ./ repmat(mags(vertind),1,3);
%                     newvertcoord = newvertcoord - (thisvertvec);
%                     thisothervec = vecsL(otherind,:) .* j ./ repmat(mags(otherind),1,3);
%                     newothercoord = newothercoord - (thisothervec);
%                 end
%             end
%             interiorL(vertind,:) = newvertcoord;
%             interiorL(otherind,:) = newothercoord;
%         end
                
        outL = Lwhite;
        outL.vertices = interiorL;
        save(outL,Linteriorsurf);
        
            
        
        
        
        
        
        vecsR = Rwhite.vertices - Rpial.vertices;
        mags = sqrt(sum((vecsR.^2),2));
        normvecsR = vecsR .* distance ./ repmat(mags,1,3);
        normvecsR(isnan(normvecsR)) = 0;
        interiorR = Rwhite.vertices + normvecsR;
        
        
        
%         white_to_interior_distances = pdist2(Rwhite.vertices,interiorR);
%         toownvert_distances = diag(white_to_interior_distances,0);
%         distance_diffs = repmat(toownvert_distances,1,size(white_to_interior_distances,2)) - white_to_interior_distances;
%         distance_diffs(mags==0,:) = 0;
%         distance_diffs(:,mags==0) = 0;
%         for i = 1:size(distance_diffs,1)
%             theseneighs = surfneighbors(i,:); theseneighs(isnan(theseneighs)) = [];
%             distance_diffs(i,theseneighs) = 0;
%         end
%         smaller_than_toownvert_inds = find(distance_diffs>.01);
%         [~, sortindex] = sort(distance_diffs(smaller_than_toownvert_inds),'ascend');
%         
%         clear distance_diffs white_to_interior_distances
%         
%         for i = 1:length(sortindex)
%             if mod(i,1000)==0
%                 disp([num2str(i) ' of ' num2str(length(sortindex))])
%             end
%             thisind = smaller_than_toownvert_inds(sortindex(i));
%             [vertind, otherind] = ind2sub([length(mags) length(mags)],thisind);
%             whitecoord = Rwhite.vertices(vertind,:);
%             newvertcoord = interiorR(vertind,:);
%             newothercoord = interiorR(otherind,:);
%             
%             for j = [.001 : .001 : 20]
%                 if pdist2(whitecoord,newothercoord) >= pdist2(whitecoord,newvertcoord)
%                     break
%                 else
%                     thisvertvec = vecsR(vertind,:) .* j ./ repmat(mags(vertind),1,3);
%                     newvertcoord = newvertcoord - (thisvertvec);
%                     thisothervec = vecsR(otherind,:) .* j ./ repmat(mags(otherind),1,3);
%                     newothercoord = newothercoord - (thisothervec);
%                 end
%             end
%             interiorR(vertind,:) = newvertcoord;
%             interiorR(otherind,:) = newothercoord;
%         end
                
        
        
        outR = Rwhite;
        outR.vertices = interiorR;
        save(outR,Rinteriorsurf);
        
    %end
    
    medial_mask_L = '/home/data/scripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii';
    medial_mask_R = '/home/data/scripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';
    
        system(['wb_command -volume-to-surface-mapping /home/data/subjects/' subject '/DTI/DTI_avg_ec_FA_MNI.nii.gz /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.midthickness.32k_fs_LR.surf.gii /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_L.func.gii -ribbon-constrained ' Linteriorsurf ' /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.L.white.32k_fs_LR.surf.gii']);
        system(['wb_command -volume-to-surface-mapping /home/data/subjects/' subject '/DTI/DTI_avg_ec_FA_MNI.nii.gz /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.midthickness.32k_fs_LR.surf.gii /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_R.func.gii -ribbon-constrained ' Rinteriorsurf ' /home/data/subjects/' subject '/fs_LR/MNI/fsaverage_LR32k/' subject '.R.white.32k_fs_LR.surf.gii']);
        system(['wb_command -cifti-create-dense-timeseries /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_LR.dtseries.nii -left-metric /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_L.func.gii -roi-left ' medial_mask_L ' -right-metric /home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_R.func.gii -roi-right ' medial_mask_R]);
        if d==1
            cifti_out = ft_read_cifti_mod(['/home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_LR.dtseries.nii']);
            cifti_out.dimord = 'scalar_pos';
        end
        thisdistdata = ft_read_cifti_mod(['/home/data/subjects/' subject '/DTI/' subject '.FA_within_' num2str(distance) 'mm_LR.dtseries.nii']);
        cifti_out.data(:,d) = thisdistdata.data;
        cifti_out.mapname{d} = ['FA within ' num2str(distance) 'mm of white/gray boundary'];
    
end

delete(['/home/data/subjects/' subject '/DTI/' subject '.FA_within_*'])
ft_write_cifti_mod(['/home/data/subjects/' subject '/DTI/FA_near_cortex_LR.dscalar.nii'],cifti_out)
