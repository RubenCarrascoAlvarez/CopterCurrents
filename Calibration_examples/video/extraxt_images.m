%
% This script extract frames from Drone videos for Camera calibration
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

dt = 0.5; % time distace in seconds between frame extracted.


vr = VideoReader('DJI_Chessboard_video_01.MP4'); % video with chessboard pattern
n = 1;
for t = 0:dt:vr.Duration
    vr.CurrentTime = t;
    img = rgb2gray(vr.readFrame());
    image_name = ['Image' num2str(n) '.tif'];
    imwrite(img,image_name);
    n = n + 1;
    
end

vr = VideoReader('DJI_Chessboard_video_02.MP4'); % video with chessboard pattern
for t = 0:dt:vr.Duration
    vr.CurrentTime = t;
    img = rgb2gray(vr.readFrame());
    image_name = ['Image' num2str(n) '.tif'];
    imwrite(img,image_name);
    n = n + 1;
    
end

% note after extract the all pictures, remove manualy the pictures where
% the chessboard pattern is not displayed completelly.
