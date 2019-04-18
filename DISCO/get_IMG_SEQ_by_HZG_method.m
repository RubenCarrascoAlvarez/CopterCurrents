function IMG_SEQ = get_IMG_SEQ_by_HZG_method(Georeference_Struct_config,heading,pitch,roll,timestamp_video,Height_drone)
% Retrieve IMG_SEQ structure using the FOV calibration described in the
% following publication:
%
% M. Streßer, R. Carrasco and J. Horstmann, "Video-Based Estimation of 
% Surface Currents Using a Low-Cost Quadcopter," in IEEE Geoscience and 
% Remote Sensing Letters, vol. 14, no. 11, pp. 2027-2031, Nov. 2017.
% doi: 10.1109/LGRS.2017.2749120
%
% The FOV projection is retrieved by gridDJIFrame_nadir.m 
%
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


% get altitude
altitude = Height_drone + Georeference_Struct_config.offset_home2water_Z;

% crete video object
v = VideoReader(Georeference_Struct_config.video_fname);

% get frames
IMG_SEQ = [];

for i1 = 1:length(Georeference_Struct_config.video_ts)

    disp(['image ' num2str(i1) ' of ' num2str(length(Georeference_Struct_config.video_ts))])
    
    % get frame from video
    v.CurrentTime = Georeference_Struct_config.video_ts(i1);
    img_RGB = readFrame(v);
    img_gray = rgb2gray(img_RGB);

    if isempty(IMG_SEQ)
        % get grid
        [X_eq,Y_eq] = gridDJIFrame_nadir(img_gray,altitude,Georeference_Struct_config.DISCO_CamCalib);
        % note: [X_eq,Y_eq] must be permuted to follow the IMG_SEQ
        % dimension requerimet: monotonic, ascending and (Nx,Ny,Nt)
        

        % create structure
        % IMG_SEQ.IMG = nan(size(X_eq,2),size(X_eq,1),length(Georeference_Struct_config.video_ts));
        IMG_SEQ.IMG = zeros(size(X_eq,2),size(X_eq,1),length(Georeference_Struct_config.video_ts),class(img_gray));
        IMG_SEQ.gridX = permute(X_eq,[2 1]);
        IMG_SEQ.gridY = permute(Y_eq,[2 1]);
        IMG_SEQ.dt = Georeference_Struct_config.dt;
        IMG_SEQ.dx = IMG_SEQ.gridX(2,1) - IMG_SEQ.gridX(1,1);
        IMG_SEQ.dy = IMG_SEQ.gridY(1,2) - IMG_SEQ.gridY(1,1);
        IMG_SEQ.ts_video = timestamp_video + (Georeference_Struct_config.video_ts/86400);
        IMG_SEQ.altitude = altitude;
        IMG_SEQ.pitch = pitch;
        IMG_SEQ.roll = roll;
        IMG_SEQ.heading = heading;
        IMG_SEQ.Georeference_Struct_config = Georeference_Struct_config;
    end

    % save image
    % IMG_SEQ.IMG(:,:,i1) = permute(img_gray,[2 1]);
    %
    % save image
    % flip Y axis + permute dim1 and dim2
    IMG_SEQ.IMG(:,:,i1) = permute(flipud(img_gray),[2 1]);

end

% destroy video Reader object
delete(v);
clear v;

end

