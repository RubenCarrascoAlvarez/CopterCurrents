%
% Addjust the FOV Calibration from a known distance.
% The Camera Calibration structure 'DISCO_CamCalib' can be used later by 
% get_IMG_SEQ_by_HZG_method.m
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

video_fname = '/media/d1/Drone_current_fit/data/Calibration/phantom3/BH/DJI_15m.MP4';

real_distance = [];


time_frame_sec = 0;
Height_offset = 0;
camera_offset_Z = 0;
DISCO_calibration_filename = []; % add full filename where the calibaration
                                 %will be stored.

disp('Measure a known distance');

[calibrated_FOV_diag,size_frame,DISCO_CamCalib] = ...
   calibration_FOV(video_fname,real_distance,time_frame_sec,Height_offset,...
   camera_offset_Z,DISCO_calibration_filename);
 
% display parameters
disp(DISCO_CamCalib)

% save calibration to file
calib_name = 'FOV_manual';
path2DISCO_Calibration_files = '/media/d1/Drone_current_fit/code/DISC0_camera_calibration/DISCO_Calibration_files/';
filename = [path2DISCO_Calibration_files calib_name '_' ...
             num2str(DISCO_CamCalib.size_X) 'x'  ...
             num2str(DISCO_CamCalib.size_Y)  '.mat'];
         
save(filename,'DISCO_CamCalib');