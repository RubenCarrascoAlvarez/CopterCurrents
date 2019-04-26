%
% Run_DISCO_script.m
%
% Run DISCO performing the current fit from a video in 5 steps:
%
% 1. Record a nadir video of the water surface by a quad-copter with a camera
% gimbal. The camera must point perpendicular to the water (nadir) and the
% Drone must be static. A 30 seconds video is enough. The Drone altitude
% is stored in the video metadata.
%
% 2. Video rectification: The images from the video are rectified to in real world
% coordinates.
%
% 3. Slice the rectified video in several squared regions, and set the fit parame-
% ters.
%
% 4. Fit the most probable current, applying the wave dispersion relation (a
% spectral energy-based maximization technique) in every squared region.
%
% 5. Filter obtained currents with Signal to Noise Ratio, and plot currents
% in Camera reference grid or UTM coordinates.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read Drone video - step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% video to be analysed
% video_fname = '/path2thevideo/Drone_video.MP4';
video_fname = '/media/d1/Drone_current_fit/data/over_the_river/PHANTOM3/20170404/DJI_0002.MP4';

% time stamps to be used in the video [initial_time  end_time]
time_limits = [5 35];

% time between frames in seconds
dt = 0.10;

% distance between the home point and the water surface in meters
offset_home2water_Z = 0.8;

% FOV calibration structure to be used
% DISCO_calibration_filename = '/path2CameraCalibration_folder/Phantom3_v1_FOV_3840x2160.mat';
DISCO_calibration_filename = '/home/carrasco/Matlab/code/DISCO/DISCO_CameraCalib/Phantom3_v1_FOV_3840x2160.mat';
% DISCO_calibration_filename = '/path2CameraCalibration_folder/Phantom3_v1_OpenCV_3840x2160.mat';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rectify video and obtain REAL WORLD coordenates - step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

% Create a Georeference configuration structure
[Georeference_Struct_config] = create_georeference_struct(...
                video_fname, dt, time_limits, offset_home2water_Z, ...
                DISCO_calibration_filename);

% Apply Georeference_Struct_config and retrieve the Image sequence
% corrected
IMG_SEQ = run_Georeference_Struct_config(Georeference_Struct_config);

% save('IMG_SEQ.mat','IMG_SEQ','-v7.3','-nocompression');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice Rectified sequence and define fit parameters - step 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% generate fit structure

% square fit size in meters
sq_size_m = []; 

% distance between square fit in meters
sq_dist_m = []; 

% 2D mask inidication the valid points for the fit
mask_2D = []; 

% percentage of area to nan to do not use the square to fit
nan_percentage_thr = 5;

% Ux current range to fit [m/s]
Ux_limits_FG = [-2.0 2.0];

% Uy current range to fit [m/s]
Uy_limits_FG = [-2.0 2.0];

% First guess step [m/s]
U_FG_res = 0.1;

% First guess filter width [rad/s]
w_width_FG = 1;

% Second guess step [m/s]
U_SG_res = U_FG_res/10;

% Second guess filter width [rad/s]
w_width_SG = w_width_FG/2;

% wavelength limits to be used in the fit 
% [shorter_waveLength longer_waveLength] in meters
waveLength_limits_m = []; 

% wave Period limits to be used in the fit 
% [smaller_wavePeriod longer_wavePeriod] in seconds
wavePeriod_limits_sec = []; 

% water depth mask in meters
% water depth in meters of every pixels in 
% water_depth_mask_2D(IMG_SEQ.gridX,IMG_SEQ.gridY)
%
% Note: if water_depth_mask_2D is a scalar, the same depth will be applied
% for the complete grid.
water_depth_mask_2D = 10; % 10 meters


% generate structure with the fit structure
STCFIT = generate_STCFIT_from_IMG_SEQ(IMG_SEQ, sq_size_m, sq_dist_m,mask_2D,nan_percentage_thr,water_depth_mask_2D,...
         Ux_limits_FG,Uy_limits_FG,U_FG_res,w_width_FG,U_SG_res,w_width_SG,waveLength_limits_m,...
         wavePeriod_limits_sec);

     
% plot STCFIT structure
[h] = plot_STCFIT(STCFIT);   
saveas(h,'STCFIT_squares_distribution.png')
close(h);

% % display  spectrum in window 'n_window'
% n_window = 38;
% display_fit_guess(IMG_SEQ,STCFIT,n_window)     


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run current fit  - step 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% run current fit in every square
STCFIT = run_current_fit(IMG_SEQ,STCFIT);

% save data
save('paper_FOV.mat','STCFIT','IMG_SEQ','-v7.3','-nocompression');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get current maps and plot results - step 5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% filter with SNR
% choose going to direction for currents
currentdir_flag = 1;
% set Signal to Noise Ratio threshold 
SNR_thr = 0; 
% set Signal to Noise Ratio density threshold 
SNR_density_thr = 6; 

% retrieve current maps
[UTM_currents, Camera_currents] = get_currents_from_STCFIT(STCFIT,SNR_thr,SNR_density_thr,currentdir_flag);

% plot in camera coordenates
% scale factor for arrows
arrow_scale = 20; 

% plot current maps in camera coordenate system
h1  = plot_currents_map(Camera_currents,STCFIT,arrow_scale);
saveas(h1,'Current_map_camera_grid.png')
close(h1);

% plot in UTM coordenates (deg2utm library needed)
if ~isempty(UTM_currents)
    h2  = plot_currents_map(UTM_currents,STCFIT,arrow_scale);
    saveas(h2,'Current_map_UTM_grid.png');
    close(h2);
end
