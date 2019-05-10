function [Georeference_Struct_config] = create_georeference_struct(...
        video_fname, dt, time_limits, offset_home2water_Z, ...
        DISCO_calibration_filename)
%Create a Georeference_Struct_config. The Georeference_Struct_config
%   defines the video parameters for the georeferencing method. 
% 
%   Input:
%     video_fname: video filename (complete path);
%     dt: distance in seconds between frames to georeference. [0.16]
%     time_limits: [initial_time_video_sec end_time_video_sec]. Vector
%                  defining the video data to be used. [0 30]
%     offset_home2water_Z: distance between home position and water surface 
%                          in meters. [0]
%     DISCO_calibration_filename: filename with the calibration parameters.
%
%
%   Output:
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

if exist('video_fname','var') == 0 || isempty(video_fname)
    video_fname = [];
end

if exist('dt','var') == 0 || isempty(dt)
    dt = 0.16;
end

if exist('time_limits','var') == 0 || isempty(time_limits)
     time_limits = [0 30];
end

if exist('offset_home2water_Z','var') == 0 || isempty(offset_home2water_Z)
    offset_home2water_Z = 0;
end

% check dt => must be multiple of 1/v.FrameRate
if ~isempty(video_fname)
    v = VideoReader(video_fname);
    dt_gap_video =  1/v.FrameRate;
    
    if rem(dt,dt_gap_video) ~= 0
        closes_dt = round(dt/dt_gap_video)*dt_gap_video;
        warning(['create_georeference_struct: dt = ' num2str(dt) ...
                 '  exchanged by: ' num2str(closes_dt)]);
        dt = closes_dt;
    end
end

% check time_limits(1)
if time_limits(1) < 0
    time_limits(1) = 0;
    warning('create_georeference_struct: time_limits(1) < 0');
    disp('time_limits(1) set to 0 seconds');
end

% check time_limits(2)
% max_duration_limit = v.Duration;
max_duration_limit = v.Duration - (1/v.FrameRate);
% max_duration_limit = v.Duration - (2/(v.FrameRate));
if time_limits(2) > max_duration_limit
    time_limits(2) = max_duration_limit;
    warning('create_georeference_struct: time_limits(2) > Video duration');
    disp(['time_limits(2) set to ' num2str(max_duration_limit) ' seconds.']);
end

% load calibration data
if isfile(DISCO_calibration_filename)
    load(DISCO_calibration_filename,'DISCO_CamCalib');
    if exist('DISCO_CamCalib','var') == 0
        error('create_georeference_struct: DISCO_calibration_filename data not valid.'); 
    end
    
else
    disp(DISCO_calibration_filename)
    error('create_georeference_struct: DISCO_calibration_filename not found.'); 
end

% check frame size

% crete video object
v = VideoReader(video_fname);

if isfield(DISCO_CamCalib,'ny') == 1
    ny = DISCO_CamCalib.ny;
else
    ny = DISCO_CamCalib.size_Y;
end

if isfield(DISCO_CamCalib,'nx') == 1
    nx = DISCO_CamCalib.nx;
else
    nx = DISCO_CamCalib.size_X;
end

% check image size
if v.Height ~= ny
    disp(['DISCO_CamCalib Height: ' num2str(ny)]);
    disp(['video Height: ' num2str(v.Height)]);
    error('video Height => Y axis (vertical) does not match with DISCO_CamCalib');
end

if v.Width ~= nx
    disp(['DISCO_CamCalib Width: ' num2str(nx)]);
    disp(['video Width: ' num2str(v.Width)]);
    error('video Width => X axis (horizontal) does not match with DISCO_CamCalib');
end


% save data to struct
Georeference_Struct_config.video_fname = video_fname;
Georeference_Struct_config.dt = dt;
Georeference_Struct_config.time_limits = time_limits;
Georeference_Struct_config.video_ts =  time_limits(1):dt:time_limits(2);
Georeference_Struct_config.offset_home2water_Z = offset_home2water_Z;
Georeference_Struct_config.DISCO_CamCalib = DISCO_CamCalib;


end

