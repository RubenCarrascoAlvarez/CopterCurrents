%
% This script extracts frames from Drone videos. This is a previous step to
% run the Camera Calibration process (Caltech toolbox calibration or OpenCV 
% calibration).
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

dt = 0.5; % gap time between frames in seconds.

n = 1; % frame name counter

% extract frames  from video 1
vr = VideoReader('DJI_Chessboard_video_01.MP4'); % video with chessboard pattern

for t = 0:dt:vr.Duration
    
    vr.CurrentTime = t; % set time in video reader
    img = rgb2gray(vr.readFrame()); % extract frame and convert to a grey scales image
    image_name = ['Image' num2str(n) '.tif']; % generate image filename
    imwrite(img,image_name); % save image in disk
    n = n + 1; % increase counter
end

% extract frames  from video 2
vr = VideoReader('DJI_Chessboard_video_02.MP4'); % video with chessboard pattern
for t = 0:dt:vr.Duration
    
    vr.CurrentTime = t; % set time in video reader
    img = rgb2gray(vr.readFrame()); % extract frame and convert to a grey scales image
    image_name = ['Image' num2str(n) '.tif']; % generate image filename
    imwrite(img,image_name); % save image in disk
    n = n + 1; % increase counter
end

% Note: after extract the all pictures, remove manually the pictures where
% the chessboard pattern is not displayed completelly. Then the pictures
% are ready to be used by the Caltech toolbox calibration or the OpenCV 
% python script.
