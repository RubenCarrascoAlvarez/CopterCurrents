%
% Addjust the FOV Calibration from a known distance.
% The Camera Calibration structure 'CopterCurrents_CamCalib' can be used later by 
% get_IMG_SEQ_by_HZG_method.m
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

video_fname = '/media/d1/Drone_current_fit/data/Calibration/phantom3/BH/DJI_15m.MP4';

% get video position in DJI drone
[CamPos_ST] = get_Camera_Position_Struct(video_fname);

% % Add data manually if is not a DJI drone 
% CamPos_ST = struct('LONGITDE',LON,'LATITUDE',LAT,'Height',Height,...
%                    'timestamp',timestamp,'yaw',yaw,'pitch',pitch,...
%                    'roll',roll,'extra',mediainfo_string);


real_distance = []; % empty value => will be asked in the command line.

time_frame_sec = 0;
Height_offset = 0;
camera_offset_Z = 0;

disp('Measure a known distance');

[calibrated_FOV_diag,size_frame,CopterCurrents_CamCalib] = ...
   calibration_FOV(video_fname,real_distance,time_frame_sec,Height_offset,...
   camera_offset_Z,[],CamPos_ST.Height);
 
% display parameters
disp(CopterCurrents_CamCalib)

% save calibration to file
  % customize your own name and CopterCurrents_CameraCalib path
calib_name = 'FOV_manual';
path2CopterCurrents_Calibration_files = '/media/d1/Drone_current_fit/githup/CopterCurrents/CopterCurrents/CopterCurrents_CameraCalib/';

filename = [path2CopterCurrents_Calibration_files calib_name '_' ...
             num2str(CopterCurrents_CamCalib.size_X) 'x'  ...
             num2str(CopterCurrents_CamCalib.size_Y)  '.mat'];
         
save(filename,'CopterCurrents_CamCalib');
disp(['Saved CameraCalib file: ' filename]);

