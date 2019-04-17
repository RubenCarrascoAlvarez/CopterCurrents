% This script process the video 'DJI_0002.MP4' used in the paper: 
% 
%   M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
%   of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
%   and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
%   doi: 10.1109/LGRS.2017.2749120
%
% The Georeferencing method use the Caltech library (compatible with OpenCV).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de>
%
% This file is part of DISCO
%
% DISCO has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
% DISCO is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% DISCO is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with DISCO.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
close all;

video_fname = '/media/d1/Drone_current_fit/data/over_the_river/PHANTOM3/20170404/DJI_0002.MP4';
dt = 0.10;
time_limits = [0 60];
offset_home2water_Z = 0;


DISCO_calibration_filename = '/media/d1/Drone_current_fit/code/DISC0_camera_calibration/DISCO_Calibration_files/Phantom3_v1_OpenCV_3840x2160.mat';

% get georeference configuration
[Georeference_Struct_config] = create_georeference_struct(...
                video_fname, dt, time_limits, offset_home2water_Z, ...
                DISCO_calibration_filename);

% geo reference data
IMG_SEQ = run_Georeference_Struct_config(Georeference_Struct_config);

% save('IMG_SEQ_0.mat','IMG_SEQ','-v7.3','-nocompression');

% generate fit structure
sq_size_m = 8; % square fit size in meters
sq_dist_m = sq_size_m/2; % distance between square fit in meters
mask_2D = []; % 2D mask inidication the valid points for the fit
nan_percentage_thr = 5; % percentage of area to nan to do not use the square to fit

% fit paramters
Ux_limits_FG = [-2.0 2.0];
Uy_limits_FG = [-2.0 2.0];
U_FG_res = 0.1;
w_width_FG = 1;
U_SG_res = U_FG_res/10;
w_width_SG = w_width_FG/2;
waveLength_limits_m = 2*pi./[10.6657 1.5641]; % [min_waveLength max_waveLength] used in the fit % K [10.6657 1.5641 ]
wavePeriod_limits_sec = 2*pi./[17 5]; % [min_wavePeriod max_wavePeriod] used in the fit             % W [5 17]
water_depth_mask_2D = [];

STCFIT = generate_STCFIT_from_IMG_SEQ(IMG_SEQ, sq_size_m, sq_dist_m,mask_2D,nan_percentage_thr,water_depth_mask_2D,...
         Ux_limits_FG,Uy_limits_FG,U_FG_res,w_width_FG,U_SG_res,w_width_SG,waveLength_limits_m,...
         wavePeriod_limits_sec);

% plot STCFIT
% [h] = plot_STCFIT(STCFIT);     

% % display paramters in window  
% n_window = 3168;
% display_fit_guess(IMG_SEQ,STCFIT,n_window)     



% run current fit in every square
STCFIT = run_current_fit(IMG_SEQ,STCFIT);


% plot STCFIT
% [h] = plot_STCFIT(STCFIT);

save('paper_caltech.mat','STCFIT','IMG_SEQ','-v7.3','-nocompression');


currentdir_flag = 1;
SNR_thr = 0;
SNR_density_thr = 3;
[UTM_currents, Camera_currents] = get_currents_from_STCFIT(STCFIT,SNR_thr,SNR_density_thr,currentdir_flag);

arrow_scale = 5;
h  = plot_currents_map(Camera_currents,STCFIT,arrow_scale);
h  = plot_currents_map(UTM_currents,STCFIT,arrow_scale);

