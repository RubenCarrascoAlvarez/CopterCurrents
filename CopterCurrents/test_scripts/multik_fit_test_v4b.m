%
% This script process the video '20170404_over_Elbe.MP4' used in the paper: 
% 
%   M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
%   of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
%   and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
%   doi: 10.1109/LGRS.2017.2749120
%
% The Georeferencing method is the FOV calibration (as in the paper).
%
%
% The video is available in: 
%   https://zenodo.org/record/3088437
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
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation 
% of Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience 
% and Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
%
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read Drone video - step 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% video to be analysed
% video_fname = '/media/d1/Miami/CC/4d_0002_shear_ADCP/DJI_0002.MP4';
video_fname = '/media/4TB_NTFS_disk1/CopterCurrents-master/videos/4b/DJI_0013.MP4'; % wind and currents in same direction

% get video position in DJI drone
[CamPos_ST] = get_Camera_Position_Struct(video_fname);

% % Add data manually if is not a DJI drone 
% CamPos_ST = struct('LONGITUDE',LON,'LATITUDE',LAT,'Height',Height,...
%                    'timestamp',timestamp,'yaw',yaw,'pitch',pitch,...
%                    'roll',roll,'extra',mediainfo_string);


% time stamps to be used in the video [initial_time  end_time]
% time_limits = [5 10];
time_limits = [1 39];

% time between frames in seconds
dt = 0.12;
% dt = 0.06;

% distance between the home point and the water surface in meters
offset_home2water_Z = 0.1;

% FOV calibration structure to be used
% CopterCurrents_calibration_filename = '/media/4TB_NTFS_disk1/CopterCurrents-master/CopterCurrents/CopterCurrents_CameraCalib/Matrice210_v1_FOV_3840x2160.mat';
CopterCurrents_calibration_filename = '/media/4TB_NTFS_disk1/CopterCurrents-master/CopterCurrents/CopterCurrents_CameraCalib/Matrice210_v1_Caltech_3840x2160.mat';

% Create a Georeference configuration structure
[Georeference_Struct_config] = create_georeference_struct(...
                video_fname, dt, time_limits, offset_home2water_Z, ...
                CopterCurrents_calibration_filename,CamPos_ST);
            
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rectify video and obtain REAL WORLD coordenates - step 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          

% Apply Georeference_Struct_config and retrieve the Image sequence
% corrected
IMG_SEQ = run_Georeference_Struct_config(Georeference_Struct_config);

% save('IMG_SEQ.mat','IMG_SEQ','-v7.3','-nocompression');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Slice Rectified sequence and define fit parameters - step 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% generate fit structure

% square fit size in meters
sq_size_m = 20; 

% distance between square fit in meters
sq_dist_m = sq_size_m; 

% 2D mask inidication the valid points for the fit
mask_2D = []; 

% percentage of area to nan to do not use the square to fit
nan_percentage_thr = 5;

% Ux current range to fit [m/s]
Ux_limits_FG = [-1.0 1.0];

% Uy current range to fit [m/s]
Uy_limits_FG = [-1.0 1.0];

% First guess step [m/s]
U_FG_res = 0.1;

% First guess filter width [rad/s]
w_width_FG = 1;

% Second guess step [m/s]
U_SG_res = 0.01;

% Second guess filter width [rad/s]
w_width_SG = w_width_FG/2;

% wavelength limits to be used in the fit 
% [shorter_waveLength longer_waveLength] in meters
% waveLength_limits_m = []; 
waveLength_limits_m = [0.15 sq_size_m/2]; 

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

% K_steps = [];
% [K_steps] = generate_K_steps(wavelength_min,wavelength_max,K_width,K_distance);
K_steps = generate_K_steps(0.15 ,2,1,0.5);
max_current_shear_expected = 0.3; % m/s


% generate structure with the fit structure
% STCFIT = generate_STCFIT_from_IMG_SEQ(IMG_SEQ, sq_size_m, sq_dist_m,mask_2D,nan_percentage_thr,water_depth_mask_2D,...
%          Ux_limits_FG,Uy_limits_FG,U_FG_res,w_width_FG,U_SG_res,w_width_SG,waveLength_limits_m,...
%          wavePeriod_limits_sec);

STCFIT_multi_K = generate_STCFIT_from_IMG_SEQ_multi_K(IMG_SEQ, sq_size_m, sq_dist_m,mask_2D,nan_percentage_thr,water_depth_mask_2D,...
         Ux_limits_FG,Uy_limits_FG,U_FG_res,w_width_FG,U_SG_res,w_width_SG,...
         waveLength_limits_m,wavePeriod_limits_sec,K_steps,max_current_shear_expected);


     
% plot STCFIT structure
[h] = plot_STCFIT(STCFIT_multi_K);   
xlabel('Distance [m]');
ylabel('Distance [m]');
set(findall(gcf,'-property','FontSize'),'FontSize',16);
print('-dpng','Window_used_n2_v4b.png');
close(h);

% % display  spectrum in window 'n_window'
% n_window = 3168;
% display_fit_guess(IMG_SEQ,STCFIT_multi_K,n_window)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run current fit and plot results - step 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% run current fit in every square
% STCFIT = run_current_fit(IMG_SEQ,STCFIT);
window_num_vec = 2;
STCFIT_multi_K = run_current_fit_multi_K(IMG_SEQ,STCFIT_multi_K,window_num_vec);

h = plot_STCFIT_multi_K(STCFIT_multi_K,window_num_vec);     

SNR_density_thr = 1.8;
wavelenth_interval = [0 inf];
polyfit_order = 2;
% fit 
close all;
% [Ux_fit,Uy_fit,depth,Ux,Uy,SNR_density,px,py] = fit_current_depth_profile(STCFIT_multi_K,window_num_vec(1),SNR_density_thr);
[Ux_fit,Uy_fit,depth,Ux,Uy,SNR_density,px,py] = fit_current_depth_profile(STCFIT_multi_K,window_num_vec(1),SNR_density_thr,wavelenth_interval,polyfit_order);
set(findall(gcf,'-property','FontSize'),'FontSize',16);
print('-dpng','Depht_fit_Wind_and_current_same_dir_v4b.png');

close all
save('multik_fit_test_v4b.mat','-v7.3','-nocompression')

% display_fit_guess(IMG_SEQ,STCFIT_multi_K,window_num_vec) 

% % save data
% save('data.mat','STCFIT','IMG_SEQ','-v7.3','-nocompression');
% 
% 
% % retrieve current maps and  plot 
% 
% % filter with SNR
% % choose going to direction for currents
% currentdir_flag = 1;
% % set Signal to Noise Ratio threshold 
% SNR_thr = 0; 
% % set Signal to Noise Ratio density threshold 
% SNR_density_thr = 3; 
% 
% % retrieve current maps
% [UTM_currents, Camera_currents] = get_currents_from_STCFIT(STCFIT,SNR_thr,SNR_density_thr,currentdir_flag);
% 
% 
% 
% % plot in camera coordenates
% % scale factor for arrows
% arrow_scale = 5; 
% 
% % plot current maps in camera coordenate system
% h  = plot_currents_map(Camera_currents,STCFIT,arrow_scale);
% 
% % plot in UTM coordenates (deg2utm library needed)
% h  = plot_currents_map(UTM_currents,STCFIT,arrow_scale);

