set name=C1-1unmixed
echo %name%
:: run curvelets first
curvelets.exe %name%.tif options_curvelets
:: run SV next
scalar_voting_2d.exe %name%_CV.tif %name%_CV_cos.mhd %name%_CV_sin.mhd options_scalar_voting_2d
:: run TV next
tensor_voting_2d.exe %name%_CV.tif %name%_CV_cos.mhd %name%_CV_sin.mhd options_scalar_voting_2d
