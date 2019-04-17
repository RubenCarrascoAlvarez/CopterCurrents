function [calibrated_FOV_diag,size_frame, DISCO_CamCalib ] = ...
          calibration_FOV(video_fname,real_distance,time_frame_sec,...
          Height_offset,camera_offset_Z,DISCO_calibration_filename)
% Make a FOV (Field of view) manual calibration.
%
%   Inputs:
%    video_fname: video file name.
%    real_distance: real distance [meters] to measure during calibration.
%                   when real_distance =[], then will be asked in the
%                   command line.
%    time_frame_sec: video time stamp to use in seconds. [0].  
%    Height_offset: Distance in meters between drone GPS antenna in home 
%                   position and the objetct to measure (or water surface). 
%                   By default [0].
%    camera_offset_Z: Distance in meters between the camera and drone GPS 
%                     antenna.
%    DISCO_calibration_filename: DISCO camera calibration filename, where
%                                the intrinsic camera parameters will be 
%                                stored. By default is empty.
%    
%    note:
%    frame_altitude = Height_GPS_drone + Height_offset + camera_offset_Z;
%
%   Outputs:
%    calibrated_FOV_diag: Diagonal FOV in degress.
%    size_frame: frame size (dim1,dim2)
%    DISCO_CamCalib: camera configuration structure
%         DISCO_CamCalib.fov_diag: diagonal field of view angle
%         DISCO_CamCalib.fov_x : x (horizontal) field of view angle;
%         DISCO_CamCalib.fov_y : y (vertical) field of view angle;
%         DISCO_CamCalib.size_X : frame pixel size in x (horizontal);
%         DISCO_CamCalib.size_Y : frame pixel size in y (vertical);
%         DISCO_CamCalib.camera_offset_Z = vertical (Z) distance in meter
%                                          between GPS and camera;
%         DISCO_CamCalib.source = 'Manual_FOV_calibration';
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


if exist('time_frame_sec','var') == 0
   time_frame_sec = 0; 
end

if exist('Height_offset','var') == 0
   Height_offset = 0; 
end

if exist('camera_offset_Z','var') == 0 || isempty(camera_offset_Z)
   camera_offset_Z = 0; 
end

if exist('DISCO_calibration_filename','var') == 0 
   DISCO_calibration_filename = []; 
end


% get metadata from video
[LON_drone,LAT_drone,Height_GPS_drone,timestamp_video, mediainfo_string,heading,pitch,roll] = get_mediainfo_information(video_fname);

% initialization vidoe object
v = VideoReader(video_fname);
v.CurrentTime = time_frame_sec;
img_RGB = readFrame(v);
img_gray = rgb2gray(img_RGB);

% note: readFrame is reading the frame:
% dim1 = vertical(Y) 
% dime = horizontal(X) 
% dim1 is fliped.
img_gray_flip = flipud(img_gray);

size_X = size(img_gray,2);
size_Y = size(img_gray,1);

% get pixel distance
x_vec_pix = 1:size_X;
y_vec_pix = 1:size_Y;
[X_pix,Y_pix] = meshgrid(x_vec_pix,y_vec_pix);

maximum_n_measures = 1;

[distance, distance_x,distance_y, initial_point,end_point ] = measure_distance_interactive(X_pix,Y_pix,img_gray_flip,maximum_n_measures);
initial_point = round(initial_point);
end_point = round(end_point);

% check real distance value
if isempty(real_distance)
    
    % ask value in command line
    prompt = 'Write distance in meters: ';
    real_distance = input(prompt);
    
    disp(['distance used: ' num2str(real_distance) ' [meters]']);
end

% get altitude
% Height_offset = 0;
altitude = Height_GPS_drone + Height_offset + camera_offset_Z;

% real distance measured
% real_distance = 10; % 10 meters (triangle side)
% real_distance = 10 * sqrt(2); % diagonal

% test FOV 
fov_diag_interval = 70:0.1:90;
distance_error = nan(size(fov_diag_interval));

% display FOV range to be evaluated
disp(['testing FOV from ' num2str(fov_diag_interval(1)) ' to ' num2str(fov_diag_interval(end))]);


for i1 = 1:length(fov_diag_interval)
    
    if rem(fov_diag_interval(i1),1) == 0
       disp(['FOV: ' num2str(fov_diag_interval(i1))]);
    end

    fov_diag = fov_diag_interval(i1);
    fov_x =  2 * atand(size_X * tand(fov_diag/2) / sqrt(size_X^2 + size_Y^2));
    fov_y = 2 * atand(size_Y * tand(fov_x/2) / size_X);

    x_vec =  linspace(altitude*tand(-fov_x/2),altitude * tand(fov_x/2),size_X);
    y_vec =  linspace(altitude*tand(-fov_y/2),altitude * tand(fov_y/2),size_Y);

    [X,Y] = meshgrid(x_vec,y_vec);

    % get error
    x_point_initial = X(initial_point(2),initial_point(1));
    y_point_initial = Y(initial_point(2),initial_point(1));

    x_point_end = X(end_point(2),end_point(1));
    y_point_end = Y(end_point(2),end_point(1));

    distance_fov =  sqrt((x_point_end-x_point_initial).^2 + (y_point_end-y_point_initial).^2);
    distance_error(i1) =  distance_fov - real_distance;

end

% interpolate best fit (Error = 0)
fov_diag_best = interp1(distance_error,fov_diag_interval,0);

figure;
plot(fov_diag_interval,distance_error,'.r'); hold on;
plot(fov_diag_best,0,'+k')
grid on;
title(['Best diagonal FOV: ' num2str(fov_diag_best)]);

% store value to return
calibrated_FOV_diag =  fov_diag_best;
size_frame = size(img_gray);

% create camera_fov_pararam structure                           
size_X = size_frame(2);
size_Y = size_frame(1);

DISCO_CamCalib.fov_diag = calibrated_FOV_diag;
DISCO_CamCalib.fov_x =  2 * atand(size_X * tand(calibrated_FOV_diag/2) / sqrt(size_X^2 + size_Y^2));
DISCO_CamCalib.fov_y =  2 * atand(size_Y * tand(calibrated_FOV_diag/2) / sqrt(size_X^2 + size_Y^2));
DISCO_CamCalib.size_X = size_X;
DISCO_CamCalib.size_Y = size_Y;
DISCO_CamCalib.camera_offset_Z = 0.0;
DISCO_CamCalib.source = 'Manual_FOV_calibration';


% save calibration if the DISCO_calibration_filename parameter is provided.
if ~isempty(DISCO_calibration_filename)
    save(DISCO_calibration_filename,'DISCO_CamCalib');
    disp(['Calibration saved: ' DISCO_calibration_filename]);
end


% display real size
[X_eq,Y_eq,IMG_eq] = gridDJIFrame_nadir(img_gray,altitude,DISCO_CamCalib);

% check measured distance
[distance2, distance_x2,distance_y2, initial_point2,end_point2 ] = measure_distance_interactive(X_eq,Y_eq,IMG_eq,inf);

% h_fig = figure('units','normalized','outerposition',[0 0 1 1]);
% h_im = pcolor(X_eq,Y_eq,IMG_eq); shading flat;
% xlabel('x distance [meters]');
% ylabel('y distance [meters]');
% h_cb = colormap(gray);
% axis xy equal tight;



end



