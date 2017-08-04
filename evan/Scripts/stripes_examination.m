function stripes_examination(filename)

a = load_untouch_nii_2D([filename '.nii.gz']);
b = demean_detrend(double(a.img));
a.img = b;
save_untouch_nii_2D(a,[filename '_dmdt.nii.gz'])
system(['mcflirt -in ' filename ' -refvol 0 -plots'])
calc_FD([filename '_mcf.par'])

