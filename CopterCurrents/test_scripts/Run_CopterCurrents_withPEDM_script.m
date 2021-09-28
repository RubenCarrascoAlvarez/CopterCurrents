%
% Run_CopterCurrents_withPEDM_script.m
%
% Run CopterCurrents performing the current fit from a video to extract the depth shear profile of the flow:
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
% 4. (optional) Fit the most probable depth uniform current approximation, applying the wave dispersion relation (a
% spectral energy-based maximization technique) in every squared region.
%
% 5. (optional) Filter obtained currents with Signal to Noise Ratio, and plot currents
% in Camera reference grid or UTM coordinates.
%
% 6. Extract the Doppler shift velocities as a function of wavenumber, and
% perform the PEDM to find the most probable shear profile.

% 7. Plot the depth profiles for a single spatial window
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (C) 2019, Ruben Carrasco <ruben.carrasco@hzg.de>
%
% This file is part of CopterCurrents
%
% CopterCurrents has been developed by Department of Radar Hydrography, at
% Helmholtz-Zentrum Geesthacht Centre for Materials and Coastal Research 
% (Germany), based in the work exposed in: 
% 
% M. StreÃŸer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
% The PEDM component used in steps 6-7 has been developed by Benjamin K.
% Smeltzer at the Norwegian University of Science and Technology
% (Trondheim, Norway). email: benjamin.smeltzer@ntnu.no

%The work is described in the following article:

%Smeltzer, B. K., Æsøy, E., Ådnøy, A.,& Ellingsen, S. Å. (2019). An improved
%method for determining near-surface currents from wave dispersion measurements. 
%Journal of Geophysical Research: Oceans, 124. https://doi.org/10.1029/2019JC015202


% CopterCurrents is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% CopterCurrents is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with CopterCurrents.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear 
close all;

addpath(genpath('C:\Sites\CopterCurrents'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read Drone video - step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% video to be analysed
%video_fname = 'G:\media\DJI_0001.mp4';
video_fname = 'G:\media\SourceVideo.MP4';
%video_fname = '/media/d1/Drone_current_fit/data/over_the_river/PHANTOM3/Zenodo/20170404_over_Elbe.MP4';

% get video position in DJI drone
[CamPos_ST] = get_Camera_Position_Struct(video_fname);

% % Add data manually if is not a DJI drone 
% CamPos_ST = struct('LONGITUDE',LON,'LATITUDE',LAT,'Height',Height,...
%                    'timestamp',timestamp,'yaw',yaw,'pitch',pitch,...
%                    'roll',roll,'extra',mediainfo_string);

% time stamps to be used in the video [initial_time  end_time]
time_limits = [1 5];

% time between frames in seconds
dt = 0.12;

% distance between the home point and the water surface in meters
offset_home2water_Z = 0.8;

% FOV calibration structure to be used
%CopterCurrents_calibration_filename = 'M:\CopterCurrents\CopterCurrents\CopterCurrents_CameraCalib/Phantom3_v1_FOV_3840x2160.mat';
%CopterCurrents_calibration_filename = '/home/carrasco/Matlab/code/CopterCurrents/CopterCurrents_CameraCalib/Phantom3_v1_FOV_3840x2160.mat';
% CopterCurrents_calibration_filename = '/path2CameraCalibration_folder/Phantom3_v1_OpenCV_3840x2160.mat';
CopterCurrents_calibration_filename = 'C:\Sites\CopterCurrents\CopterCurrents\CopterCurrents_CameraCalib/Phantom3_v1_FOV_3840x2160.mat';
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rectify video and obtain REAL WORLD coordenates - step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

% Create a Georeference configuration structure
[Georeference_Struct_config] = create_georeference_struct(...
                video_fname, dt, time_limits, offset_home2water_Z, ...
                CopterCurrents_calibration_filename,CamPos_ST);

% Apply Georeference_Struct_config and retrieve the Image sequence
% corrected
IMG_SEQ = run_Georeference_Struct_config(Georeference_Struct_config);

%%
for i = 1%:size(IMG_SEQ.IMG,3)
    %figure(1);imagesc(IMG_SEQ.IMG(:,:,i));axis image;colorbar;
    figure(2);imagesc(IMG_SEQ.gridY);axis image;colorbar;
    drawnow;
end


% save('IMG_SEQ.mat','IMG_SEQ','-v7.3','-nocompression');
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice Rectified sequence and define fit parameters - step 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% generate fit structure

% square fit size in meters
sq_size_m = [20]; 

% distance between square fit in meters
sq_dist_m = [25]; 

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
% saveas(h,'STCFIT_squares_distribution.png')
% close(h);

% % display  spectrum in window 'n_window'
% n_window = 38;
% display_fit_guess(IMG_SEQ,STCFIT,n_window)     
%% OPTIONAL (STEPS 4-5) %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run current fit  - step 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% run current fit in every square
% tic;
STCFIT = run_current_fit(IMG_SEQ,STCFIT);
% toc;


% % save data
% save('paper_FOV.mat','STCFIT','IMG_SEQ','-v7.3','-nocompression');

% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Get current maps and plot results - step 5
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% % filter with SNR
% % choose going to direction for currents
 currentdir_flag = 1;
% % set Signal to Noise Ratio threshold 
 SNR_thr = 0; 
% % set Signal to Noise Ratio density threshold 
 SNR_density_thr = 3; 
% 
% % retrieve current maps
[UTM_currents, Camera_currents] = get_currents_from_STCFIT(STCFIT,SNR_thr,SNR_density_thr,currentdir_flag);
% 
% % plot in camera coordenates
% % scale factor for arrows
 arrow_scale = 20; 
% 
% % plot current maps in camera coordenate system
h1  = plot_currents_map(Camera_currents,STCFIT,arrow_scale);
% saveas(h1,'Current_map_camera_grid.png')
% close(h1);
% 
% % plot in UTM coordenates (deg2utm library needed)
% if ~isempty(UTM_currents)
%     h2  = plot_currents_map(UTM_currents,STCFIT,arrow_scale);
%     saveas(h2,'Current_map_UTM_grid.png');
%     close(h2);
% end
% 
 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the wavenumber-dependent Doppler shift velocities and perform the PEDM - step 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define parameters specific to the k-dependent Doppler shift extraction and
%PEDM

%List of wavenumber values at which to extract Doppler shift velocities
%[rad/m]. If empty, values will be chosen between the values of
%STCFIT.fit_param.K_limits steps of the k pixel spacing of each window (assumed to be
%square)
dk = 2*pi/STCFIT.Windows.sq_size_m;% wavenumber resolution of spectrum in each spatial window
wavenumbers = dk*8:dk:10;

% (optional): whether to include 2nd harmonic of the spectrum in the fit (false by default)
include2ndHarmonic = [];

% (optional): whether to do the fit in log space (false by default)
logFlag = [];

% (optional) omegaWidthFun: function handle as a function of wavenumber i.e.
%@(k) f(k)...., specifying frequency width of the weighting function in
%frequency-angle space (constant wavenumber). Width is half-width 1/e^2
%point of a Gaussian function.
omegaWidthFun = [];

%The following OPTIONAL parameters involve post-processing of the Doppler shifts:
%SNR_filter: whether to use a signal-to-noise filter (false by default)
SNR_filter = 1.0;

%SNR_threshold: threshold signal-to-noise value for above filter (set to 2.0 by default)
SNR_threshold = 2.0;

%Peak_filter: whether to use a multiple peaks filter (false by default)
Peak_filter = 1.0;

%Peak_threshold: peak threshold of maximum value (0.5 by default)
Peak_threshold = 0.5;

%Outlier_filter: whether to use an outlier filter (quartile-based) (false by default)
Outlier_filter = 1.0;

STCFIT = generate_STCFIT_for_NSPP(STCFIT,wavenumbers,include2ndHarmonic,logFlag,...
    omegaWidthFun,SNR_filter,SNR_threshold,Peak_filter,Peak_threshold,Outlier_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now extract the Doppler shifts and run the pedm for each spatial window,
%or a subset of windows.
windowList = 6;%[5,10];%[5,10];%Subset of windows to run (set to [] or omit as input below to run over all windows).

tic;
STCFIT = run_current_fit_depth_profile(IMG_SEQ,STCFIT,windowList);
toc;


%% Step 7: plot depth profiles for an individual window.

i = 1;%Window index, indexes into 'windowList' above.

plot_depth_profile_results(STCFIT.out_depth_profile(i),3);
plot_pedm_results(STCFIT.out_depth_profile(i),4);




