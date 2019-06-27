function IMG_SEQ = get_IMG_SEQ_by_caltech_library_v2(Georeference_Struct_config,heading,pitch,roll,timestamp_video,Height_drone)
% Retrieve IMG_SEQ structure using the 'Camera Calibration Toolbox for Matlab'.
% The IMG_SEQ structure stores the video data in real size.
%   
%   Note: the caltech library 'Camera Calibration Toolbox for Matlab' is 
%   an external libary, written by Jean-Yves Bouguet,
%   jean-yves.bouguet@intel.com. This library is not part of CopterCurrents, and it
%   is available in:
%
%   http://www.vision.caltech.edu/bouguetj/calib_doc/
%   
%   The 'Camera Calibration Toolbox for Matlab' is compatible with OpenCV
%   camera calibration. This is the main reason to include this function in
%   CopterCurrents.
%  
%   If you want to retrieve the IMG_SEQ structure without any CopterCurrents 
%   external library, then use 'get_IMG_SEQ_by_HZG_method.m'.
% 
%   Input:
%     Georeference_Struct_config: Georeference configuration 
%     heading: video metadata heading (yaw) in degrees 
%     pitch: video metadata pitch in degrees 
%     roll: video metadata roll in degrees 
%     timestamp_video: time video information in seconds
%     Height_drone: distance from camera to water surface in meters.
%
%   Output:
%     IMG_SEQ structure
%     IMG_SEQ.IMG: 3D Intensty data (Nx,Ny,Nt)
%     IMG_SEQ.gridX: 2D grid X data in meters (Nx,Ny)
%     IMG_SEQ.gridY: 2D grid Y data in meters (Nx,Ny)
%     IMG_SEQ.dx: X resolution in meters.
%     IMG_SEQ.dy: Y resolution in meters.
%     IMG_SEQ.dt: time resolution in seconds.
%     IMG_SEQ.ts_video: video time stamps in datenum format. 1xNt vector.
%     IMG_SEQ.altitude: altitude to the water surface in meters.
%                       altitude  = video altitude + offset_home2water_Z
%     IMG_SEQ.pitch: video pitch (for Nadir pitch = -90)
%     IMG_SEQ.roll: video roll (for Nadir roll = 0)
%     IMG_SEQ.heading: video heading to North (yaw)
%     IMG_SEQ.Longitude: video longitude in degrees
%     IMG_SEQ.Latitude: video latitude in degrees
%     IMG_SEQ.mediainfo_string: raw mediainfo string
%     IMG_SEQ.Georeference_Struct_config: Georeference_Struct used
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


%
%     Georeference_Struct_config struct
%     Georeference_Struct_config.video_fname: video filename (complete path);
%     Georeference_Struct_config.dt: distance in seconds between frames to 
%                georeference. [initial_time_video_sec end_time_video_sec]
%     Georeference_Struct_config.time_limits: Vector defining the video 
%                                 data to be used.
%     Georeference_Struct_config.video_ts: video time stamps in detenum format
%     Georeference_Struct_config.offset_home2water_Z: distance between 
%                           home position and water surface in meters.
%     Georeference_Struct_config.CopterCurrents_CamCalib: Camera calibration
%                                                structure.
%     Georeference_Struct_config.out_resolution: resolution in meters of 
%           georeferenced data, only used when georeference_mode = 0 
%     

% get sensor offset
v = VideoReader(Georeference_Struct_config.video_fname);
% img_RGB = readFrame(v);


% get altitute
altitude = Height_drone + Georeference_Struct_config.offset_home2water_Z + ...
           Georeference_Struct_config.CopterCurrents_CamCalib.camera_offset_Z;

% get equidistant grid
v.CurrentTime = Georeference_Struct_config.video_ts(1);
img_RGB = readFrame(v);
img_gray = rgb2gray(img_RGB);

% get conversion to monotonic grid structure using Caltech library 
% [~,conv_monotonic_grid_ST] = georeferenceDJIFrame_byCaltech_EQ(img_gray,altitude, pitch,roll,Georeference_Struct_config.CopterCurrents_CamCalib);
[~,conv_monotonic_grid_ST] = georeferenceDJIFrame_byCaltech_EQ_v2(img_gray,altitude, pitch,roll,Georeference_Struct_config.CopterCurrents_CamCalib);
       
IMG_SEQ = [];

% axis_flag = 0; => georeferenceDJIFrame_byCaltech_EQ_v2
% axis_flag = 1; => georeferenceDJIFrame_byCaltech_EQ
axis_flag = 0; 
% axis flag to exchange X and Y axis, 
% and exchange dimesion 1 by 2 (IMG_SEQ format)


for i1 = 1:length(Georeference_Struct_config.video_ts)

    disp(['image ' num2str(i1) ' of ' num2str(length(Georeference_Struct_config.video_ts))])
    
    % get frame from video
    v.CurrentTime = Georeference_Struct_config.video_ts(i1);
    img_RGB = readFrame(v);
    img_gray = rgb2gray(img_RGB);
    
    % get not equidisntant grid
    [X_eq,Y_eq,IMG_eq] = apply_conv_monotonic_grid_ST(conv_monotonic_grid_ST,img_gray,axis_flag);
    % convert to single to save memory
    IMG_eq = single(IMG_eq);

    if isempty(IMG_SEQ)

        % IMG_SEQ.IMG = nan(size(X_eq,1),size(X_eq,2),length(Georeference_Struct_config.video_ts),class(IMG_eq));
        IMG_SEQ.IMG = zeros(size(X_eq,1),size(X_eq,2),length(Georeference_Struct_config.video_ts),class(img_gray));
        IMG_SEQ.gridX = X_eq;
        IMG_SEQ.gridY = Y_eq;
        IMG_SEQ.dt = Georeference_Struct_config.dt;
        IMG_SEQ.dx = conv_monotonic_grid_ST.dxdy;
        IMG_SEQ.dy = conv_monotonic_grid_ST.dxdy;
        IMG_SEQ.ts_video = timestamp_video + (Georeference_Struct_config.video_ts/86400);
        IMG_SEQ.altitude = altitude;
        IMG_SEQ.pitch = pitch;
        IMG_SEQ.roll = roll;
        IMG_SEQ.heading = heading;
        IMG_SEQ.Georeference_Struct_config = Georeference_Struct_config;
    end
    
    % save georeference picture
    IMG_SEQ.IMG(:,:,i1) = IMG_eq;

    % figure;pcolor(IMG_SEQ.gridX,IMG_SEQ.gridY,IMG_SEQ.IMG(:,:,i1)); shading flat; axis xy equal tight;

end

% destroy video Reader object
delete(v);
clear v;

end

