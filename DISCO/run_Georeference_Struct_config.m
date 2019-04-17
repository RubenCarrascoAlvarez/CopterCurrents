function [IMG_SEQ] = run_Georeference_Struct_config(Georeference_Struct_config)
% Retrieve IMG_SEQ structure from Georeference_Struct_config
% The IMG_SEQ structure contain the projected drone 3D data (x,y,time) in
% real size, according to the drone height.
%
%   Input:
%     Georeference_Struct_config struct
%     Georeference_Struct_config.video_fname: video filename (complete path);
%     Georeference_Struct_config.dt: distance in seconds between frames to 
%                georeference. [initial_time_video_sec end_time_video_sec]
%     Georeference_Struct_config.time_limits: Vector defining the video 
%                                 data to be used.
%     Georeference_Struct_config.video_ts: video time stamps in detenum format
%     Georeference_Struct_config.offset_home2water_Z: distance between 
%                           home position and water surface in meters.
%     Georeference_Struct_config.DISCO_CamCalib: Camera calibration
%                                                structure.
%     
%    
%   Output:
%     IMG_SEQ structure
%     IMG_SEQ.IMG: 3D Intensty data (Nx,Ny,Nt)
%     IMG_SEQ.gridX: 2D grid X data in meters (Nx,Ny) monotonic in ascending order.
%     IMG_SEQ.gridY: 2D grid Y data in meters (Nx,Ny) monotonic in ascending order.
%     IMG_SEQ.dx: X resolution in meters.
%     IMG_SEQ.dy: Y resolution in meters.
%     IMG_SEQ.dt: time resolution in seconds.
%     IMG_SEQ.ts_video: video time stamps in datenum format. 
%                       1xNt vector monotonic in ascending order.
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

% get metadata from video
[LON_drone,LAT_drone,Height_drone,timestamp_video, mediainfo_string,heading,pitch,roll] = get_mediainfo_information(Georeference_Struct_config.video_fname);

% check DISCO_CamCalib type
[GeoMode] = getGeoMode_from_DISCO_CamCalib(Georeference_Struct_config.DISCO_CamCalib);

% georeference data
if GeoMode == 0 % use caltech geoference
    % IMG_SEQ = get_IMG_SEQ_by_caltech_library(Georeference_Struct_config,heading,pitch,roll,timestamp_video,Height_drone); 
    IMG_SEQ = get_IMG_SEQ_by_caltech_library_v2(Georeference_Struct_config,heading,pitch,roll,timestamp_video,Height_drone); 
else  % use FOV georeference mode
    IMG_SEQ = get_IMG_SEQ_by_HZG_method(Georeference_Struct_config,heading,pitch,roll,timestamp_video,Height_drone);
end

% add lon lat and media info information
IMG_SEQ.Longitude = LON_drone;
IMG_SEQ.Latitude = LAT_drone;
IMG_SEQ.mediainfo_string = mediainfo_string;

end

