function label_from_metric(metricfile,labelcolornamefile,hem)
%function label_from_metric(metricfile,labelcolornamefile,hem)

labeltemplatefile = ['/data/cn4/laumannt/FinalLabels/Power_Neuron11_dil.' hem '.32k_fs_LR.label.gii'];

label = gifti(labeltemplatefile);

[index R G B a labelID] = textread(labelcolornamefile,'%f %f %f %f %f %s');

metric = gifti(metricfile);

label.cdata = uint32(metric.cdata);

label.private.data{1}.attributes.ArrayIndexingOrder = 'RowMajorOrder';

label.private.data{1}.attributes.Dim = max(label.private.data{1}.attributes.Dim);

label.private.metadata(5).value = ['Created from ' metricfile];

label.private.metadata(3).value = date;

label.private.label.name = labelID';

label.private.label.key = index';

label.private.label.rgba = [R G B a];

label.private.data{1}.metadata(2).value = ' ';
label.private.data{1}.metadata(1).value = ' ';

save(label,[metricfile(1:end-9) '.label.gii']);


