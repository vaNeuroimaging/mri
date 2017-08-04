
ciftifilename = 'Cluster_probability_maps_sorted_10mm_mediansizeoutlines_40sub_notintemplate_plusextended_matchedtoHCP_matchedto108_all';%.dtseries.nii';

cifti = ft_read_cifti_mod([ciftifilename '.dtseries.nii']);


networknames = {'Default','Visual','FrontoPar','FrontoPar2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MotorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'};
colors = [1 0 0;0 0 .6;1 1 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;1 1 .8;0 .4 0;.25 .25 .25];

hemnames = {'L','R'};

hemverts = 32492;

for h = 1:length(hemnames)
    
    hemdata = zeros(hemverts,1);
    heminds{h} = cifti.brainstructure((1:hemverts)+(hemverts*(h-1))) >0;
    hemdata(heminds{h}) = cifti.data((1:nnz(heminds{h})) + (nnz(heminds{1}) * (h-1)));
    
    foci = find(hemdata>0);
    
    IDs = hemdata(foci);
        
    RGB = colors(IDs,:);
    
    clear fociname classname
    for i = 1:length(foci)
        fociname{i} = ['Foci' num2str(i)];
        classname{i} = networknames{IDs(i)};
    end
    
    focifilename = [ciftifilename '_' hemnames{h} '.foci'];

    focifilemaker_workbench(foci,RGB,fociname,classname,focifilename,hemnames{h});
end