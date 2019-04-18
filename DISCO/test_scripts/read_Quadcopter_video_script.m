%
% This is a testing script to check if a Quadrocopter video can be read
% by the matlab videoreader, and if the video metadata information is  
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

close all;
clear 
clc;

video_fname = 'path to test video/test_video.MP4';
% video_fname = '/media/d1/Drone_current_fit/data/over_the_river/PHANTOM3/20170404/DJI_0002.MP4';


% initialization video object
v = VideoReader(video_fname);
v.CurrentTime = v.duration/2;
img_RGB = readFrame(v);
img_gray = rgb2gray(img_RGB);

% display frame for 5 seconds
h = figure;
colormap(gray);
imagesc(img_gray);
axis xy equal tight;
xlabel('pixels');
ylabel('pixels');
pause(5);
close(h)

% get metadata from video
[LON_drone,LAT_drone,Height_GPS_drone,timestamp_video, mediainfo_string,...
    heading,pitch,roll] = get_mediainfo_information(video_fname);

% display raw metadata info
% disp(mediainfo_string);

disp('Metadata info obtained')
disp(['Longitude :' num2str(LON_drone) ' deg']);
disp(['Latitude :' num2str(LAT_drone) ' deg']);
disp(['Altitude :' num2str(Height_GPS_drone) ' meters']);
disp(['time stamp: ' datestr(timestamp_video)]);
disp(['Camera heading:' num2str(heading) ' deg']);
disp(['Camera pitch:' num2str(pitch) ' deg']);
disp(['Camera roll :' num2str(roll) ' deg']);
