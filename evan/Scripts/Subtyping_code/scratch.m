effsizedouble = cell(3,1);
effsizesign = cell(3,1);
effsizesplitvar = cell(5,1);
effsizedouble{1} = zeros(600,1);
effsizesign{1} = zeros(600,1);
group1_subjectsize = 60;
group2_subjectsize = 60;
datasplit = 0.8;
nrepsCI = 1000;
effsizedouble{2} = zeros(600,1) + 0.5;
effsizedouble{3} = zeros(600,1) + 0.5;
effsizesign{2} = zeros(600,1) + 0.5;
effsizesign{3} = zeros(600,1) + 0.5;
effsizesplitvar{1} = zeros(600,1);
effsizesplitvar{2} = zeros(150,1) + 0.5;
effsizesplitvar{3} = zeros(150,1) + 0.5;
effsizesplitvar{4} = zeros(150,1) + 0.5;
effsizesplitvar{5} = zeros(150,1) + 0.5;
nrepsPM = 1;
[group1_data_double, group2_data_double] = SimulateTwoBiGroupDataByDoubling(group1_subjectsize,group2_subjectsize,effsizedouble{1},effsizedouble{2},effsizedouble{3});
[group1_data_sign, group2_data_sign] = SimulateTwoBiGroupDataBySign(group1_subjectsize,group2_subjectsize,effsizesign{1},effsizesign{2},effsizesign{3});
[group1_data_split, group2_data_split] = SimulateTwoBiGroupDataSplitVar(group1_subjectsize,group2_subjectsize,effsizesplitvar{1},effsizesplitvar{2},effsizesplitvar{3},effsizesplitvar{4},effsizesplitvar{5});
