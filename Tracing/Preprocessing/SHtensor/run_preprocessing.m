%run_preprocessing
farsight_trunk = 'C:\Users\arun\Research\Farsight\src\trunk';
addpath([ farsight_trunk '\Tracing\Preprocessing\SHtensor\2D']);
addpath([ farsight_trunk '\Tracing\Preprocessing\SHtensor\3D']);
addpath([ farsight_trunk '\Tracing\Preprocessing\CurveletFiltering\fdct_wrapping_matlab']);
im = readim('C:\Users\arun\Research\TensorVoting\SHtensor\Diadem\MAX_hippocampal_cropped.tif');


[out,gradmag] = scalar_voting_main(im);

opt.BlackWhite = false;
opt.FrangiScaleRange = [ 1 5];
[vesselness, scale, dir] = FrangiFilter2D(out,opt);
[vesselness2, scale, dir] = FrangiFilter2D(gradmag,opt);

figure, imagesc(vesselness);
figure, imagesc(vesselness2);