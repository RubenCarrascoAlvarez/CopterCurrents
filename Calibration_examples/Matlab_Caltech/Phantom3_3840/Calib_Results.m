% Intrinsic and Extrinsic Camera Parameters
%
% This script file can be directly executed under Matlab to recover the camera intrinsic and extrinsic parameters.
% IMPORTANT: This file contains neither the structure of the calibration objects nor the image coordinates of the calibration points.
%            All those complementary variables are saved in the complete matlab data file Calib_Results.mat.
% For more information regarding the calibration model visit http://www.vision.caltech.edu/bouguetj/calib_doc/


%-- Focal length:
fc = [ 2351.019539425834409 ; 2353.014635737576100 ];

%-- Principal point:
cc = [ 1938.717863090774244 ; 1063.300283219928133 ];

%-- Skew coefficient:
alpha_c = 0.000000000000000;

%-- Distortion coefficients:
kc = [ -0.135903056502948 ; 0.119293592777893 ; 0.000514683001335 ; 0.000996817710125 ; 0.000000000000000 ];

%-- Focal length uncertainty:
fc_error = [ 40.777158480368207 ; 40.757374791949672 ];

%-- Principal point uncertainty:
cc_error = [ 11.990496956036832 ; 8.579324031896839 ];

%-- Skew coefficient uncertainty:
alpha_c_error = 0.000000000000000;

%-- Distortion coefficients uncertainty:
kc_error = [ 0.006757131963986 ; 0.013766303406413 ; 0.000429520258092 ; 0.000670660109478 ; 0.000000000000000 ];

%-- Image size:
nx = 3840;
ny = 2160;


%-- Various other variables (may be ignored if you do not use the Matlab Calibration Toolbox):
%-- Those variables are used to control which intrinsic parameters should be optimized

n_ima = 28;						% Number of calibration images
est_fc = [ 1 ; 1 ];					% Estimation indicator of the two focal variables
est_aspect_ratio = 1;				% Estimation indicator of the aspect ratio fc(2)/fc(1)
center_optim = 1;					% Estimation indicator of the principal point
est_alpha = 0;						% Estimation indicator of the skew coefficient
est_dist = [ 1 ; 1 ; 1 ; 1 ; 0 ];	% Estimation indicator of the distortion coefficients


%-- Extrinsic parameters:
%-- The rotation (omc_kk) and the translation (Tc_kk) vectors for every calibration image and their uncertainties

%-- Image #1:
omc_1 = [ -2.296311e+00 ; -2.119041e+00 ; -4.210171e-02 ];
Tc_1  = [ -5.608740e+02 ; -3.285535e+02 ; 8.633096e+02 ];
omc_error_1 = [ 1.807685e-03 ; 1.455608e-03 ; 3.712612e-03 ];
Tc_error_1  = [ 4.519478e+00 ; 3.276285e+00 ; 1.480611e+01 ];

%-- Image #2:
omc_2 = [ -2.293592e+00 ; -2.123217e+00 ; -4.364558e-02 ];
Tc_2  = [ -5.745177e+02 ; -3.482297e+02 ; 9.837198e+02 ];
omc_error_2 = [ 2.120858e-03 ; 1.695476e-03 ; 4.208922e-03 ];
Tc_error_2  = [ 5.135305e+00 ; 3.718573e+00 ; 1.691168e+01 ];

%-- Image #3:
omc_3 = [ -2.290162e+00 ; -2.127606e+00 ; -4.147417e-02 ];
Tc_3  = [ -5.874687e+02 ; -3.600484e+02 ; 1.073508e+03 ];
omc_error_3 = [ 2.365695e-03 ; 1.883871e-03 ; 4.649406e-03 ];
Tc_error_3  = [ 5.588254e+00 ; 4.044171e+00 ; 1.849325e+01 ];

%-- Image #4:
omc_4 = [ -2.293434e+00 ; -2.120884e+00 ; -4.379817e-02 ];
Tc_4  = [ -6.161582e+02 ; -3.418646e+02 ; 1.171618e+03 ];
omc_error_4 = [ 2.654154e-03 ; 2.123584e-03 ; 5.226849e-03 ];
Tc_error_4  = [ 6.084441e+00 ; 4.405820e+00 ; 2.021846e+01 ];

%-- Image #5:
omc_5 = [ -2.329674e+00 ; -2.082796e+00 ; -4.658900e-02 ];
Tc_5  = [ -6.452006e+02 ; -2.934738e+02 ; 1.273574e+03 ];
omc_error_5 = [ 2.993420e-03 ; 2.358156e-03 ; 5.936464e-03 ];
Tc_error_5  = [ 6.603713e+00 ; 4.777006e+00 ; 2.202039e+01 ];

%-- Image #6:
omc_6 = [ -2.361145e+00 ; -2.059925e+00 ; -3.464283e-02 ];
Tc_6  = [ -6.380785e+02 ; -2.667002e+02 ; 1.263708e+03 ];
omc_error_6 = [ 2.955339e-03 ; 2.290938e-03 ; 5.862172e-03 ];
Tc_error_6  = [ 6.541020e+00 ; 4.719421e+00 ; 2.188502e+01 ];

%-- Image #7:
omc_7 = [ -2.364509e+00 ; -2.059068e+00 ; -3.085419e-02 ];
Tc_7  = [ -6.227954e+02 ; -2.328973e+02 ; 1.164019e+03 ];
omc_error_7 = [ 2.617499e-03 ; 2.079672e-03 ; 5.180901e-03 ];
Tc_error_7  = [ 6.026820e+00 ; 4.347690e+00 ; 2.014960e+01 ];

%-- Image #8:
omc_8 = [ -2.366610e+00 ; -2.063989e+00 ; -2.475968e-02 ];
Tc_8  = [ -6.120812e+02 ; -2.247755e+02 ; 1.023496e+03 ];
omc_error_8 = [ 2.210666e-03 ; 1.761870e-03 ; 4.320197e-03 ];
Tc_error_8  = [ 5.310545e+00 ; 3.826581e+00 ; 1.767064e+01 ];

%-- Image #9:
omc_9 = [ -2.363003e+00 ; -2.062633e+00 ; -2.733622e-02 ];
Tc_9  = [ -6.159053e+02 ; -2.309543e+02 ; 8.778962e+02 ];
omc_error_9 = [ 1.860975e-03 ; 1.458194e-03 ; 3.656952e-03 ];
Tc_error_9  = [ 4.574752e+00 ; 3.302604e+00 ; 1.509385e+01 ];

%-- Image #10:
omc_10 = [ NaN ; NaN ; NaN ];
Tc_10  = [ NaN ; NaN ; NaN ];
omc_error_10 = [ NaN ; NaN ; NaN ];
Tc_error_10  = [ NaN ; NaN ; NaN ];

%-- Image #11:
omc_11 = [ NaN ; NaN ; NaN ];
Tc_11  = [ NaN ; NaN ; NaN ];
omc_error_11 = [ NaN ; NaN ; NaN ];
Tc_error_11  = [ NaN ; NaN ; NaN ];

%-- Image #12:
omc_12 = [ NaN ; NaN ; NaN ];
Tc_12  = [ NaN ; NaN ; NaN ];
omc_error_12 = [ NaN ; NaN ; NaN ];
Tc_error_12  = [ NaN ; NaN ; NaN ];

%-- Image #13:
omc_13 = [ NaN ; NaN ; NaN ];
Tc_13  = [ NaN ; NaN ; NaN ];
omc_error_13 = [ NaN ; NaN ; NaN ];
Tc_error_13  = [ NaN ; NaN ; NaN ];

%-- Image #14:
omc_14 = [ NaN ; NaN ; NaN ];
Tc_14  = [ NaN ; NaN ; NaN ];
omc_error_14 = [ NaN ; NaN ; NaN ];
Tc_error_14  = [ NaN ; NaN ; NaN ];

%-- Image #15:
omc_15 = [ NaN ; NaN ; NaN ];
Tc_15  = [ NaN ; NaN ; NaN ];
omc_error_15 = [ NaN ; NaN ; NaN ];
Tc_error_15  = [ NaN ; NaN ; NaN ];

%-- Image #16:
omc_16 = [ NaN ; NaN ; NaN ];
Tc_16  = [ NaN ; NaN ; NaN ];
omc_error_16 = [ NaN ; NaN ; NaN ];
Tc_error_16  = [ NaN ; NaN ; NaN ];

%-- Image #17:
omc_17 = [ NaN ; NaN ; NaN ];
Tc_17  = [ NaN ; NaN ; NaN ];
omc_error_17 = [ NaN ; NaN ; NaN ];
Tc_error_17  = [ NaN ; NaN ; NaN ];

%-- Image #18:
omc_18 = [ NaN ; NaN ; NaN ];
Tc_18  = [ NaN ; NaN ; NaN ];
omc_error_18 = [ NaN ; NaN ; NaN ];
Tc_error_18  = [ NaN ; NaN ; NaN ];

%-- Image #19:
omc_19 = [ NaN ; NaN ; NaN ];
Tc_19  = [ NaN ; NaN ; NaN ];
omc_error_19 = [ NaN ; NaN ; NaN ];
Tc_error_19  = [ NaN ; NaN ; NaN ];

%-- Image #20:
omc_20 = [ NaN ; NaN ; NaN ];
Tc_20  = [ NaN ; NaN ; NaN ];
omc_error_20 = [ NaN ; NaN ; NaN ];
Tc_error_20  = [ NaN ; NaN ; NaN ];

%-- Image #21:
omc_21 = [ -2.379704e+00 ; -1.902355e+00 ; 4.498918e-01 ];
Tc_21  = [ -6.289579e+02 ; -1.593263e+02 ; 1.150102e+03 ];
omc_error_21 = [ 2.668546e-03 ; 1.918420e-03 ; 8.809889e-03 ];
Tc_error_21  = [ 5.879605e+00 ; 4.213984e+00 ; 1.987237e+01 ];

%-- Image #22:
omc_22 = [ -2.468724e+00 ; -1.792107e+00 ; 4.179941e-01 ];
Tc_22  = [ -6.600684e+02 ; -7.183560e+01 ; 1.257912e+03 ];
omc_error_22 = [ 2.725460e-03 ; 1.901318e-03 ; 8.787468e-03 ];
Tc_error_22  = [ 6.410779e+00 ; 4.599119e+00 ; 2.181385e+01 ];

%-- Image #23:
omc_23 = [ -2.470113e+00 ; -1.786337e+00 ; 4.036991e-01 ];
Tc_23  = [ -5.901646e+02 ; -5.205669e+01 ; 1.313126e+03 ];
omc_error_23 = [ 2.786860e-03 ; 2.052367e-03 ; 8.951661e-03 ];
Tc_error_23  = [ 6.694233e+00 ; 4.776872e+00 ; 2.281885e+01 ];

%-- Image #24:
omc_24 = [ -2.324845e+00 ; -1.961770e+00 ; 4.628683e-01 ];
Tc_24  = [ -4.319493e+02 ; -1.499050e+02 ; 1.346240e+03 ];
omc_error_24 = [ 2.879751e-03 ; 2.403235e-03 ; 9.460758e-03 ];
Tc_error_24  = [ 6.899715e+00 ; 4.892967e+00 ; 2.324167e+01 ];

%-- Image #25:
omc_25 = [ -2.170576e+00 ; -2.122350e+00 ; 5.151014e-01 ];
Tc_25  = [ -3.430948e+02 ; -2.150628e+02 ; 1.338499e+03 ];
omc_error_25 = [ 2.842364e-03 ; 2.617835e-03 ; 9.437094e-03 ];
Tc_error_25  = [ 6.865894e+00 ; 4.888981e+00 ; 2.287737e+01 ];

%-- Image #26:
omc_26 = [ -2.127508e+00 ; -2.166262e+00 ; 5.296122e-01 ];
Tc_26  = [ -3.359625e+02 ; -2.237496e+02 ; 1.326457e+03 ];
omc_error_26 = [ 2.835667e-03 ; 2.654368e-03 ; 9.425466e-03 ];
Tc_error_26  = [ 6.805549e+00 ; 4.857650e+00 ; 2.260701e+01 ];

%-- Image #27:
omc_27 = [ -2.121127e+00 ; -2.178299e+00 ; 5.403385e-01 ];
Tc_27  = [ -3.714802e+02 ; -2.248087e+02 ; 1.288314e+03 ];
omc_error_27 = [ 2.855784e-03 ; 2.578278e-03 ; 9.497609e-03 ];
Tc_error_27  = [ 6.623133e+00 ; 4.736493e+00 ; 2.193560e+01 ];

%-- Image #28:
omc_28 = [ -2.146067e+00 ; -2.155444e+00 ; 5.310482e-01 ];
Tc_28  = [ -4.419289e+02 ; -2.250593e+02 ; 1.190758e+03 ];
omc_error_28 = [ 2.786886e-03 ; 2.381173e-03 ; 9.364278e-03 ];
Tc_error_28  = [ 6.127547e+00 ; 4.386006e+00 ; 2.032789e+01 ];

