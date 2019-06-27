function [X,Y,IMG] = gridDJIFrame_nadir(img,altitude,config_str)
% Retrieve X and Y 2D coordenates using the FOV method, from DJI video.
%   Note: the video must be recorded in Nadir postion, camera perpedicular
%   to the water surface.
%     
%   Input:
%     img: grey scale raw 2D Intensity image (Ny,Nx), extrated from video  
%          according to matlab image standar:
%          dim1 => Y (vertical) data.
%          dim2 => X (Horizontal) data.
%     altitude: Camera altitude in meters to the water surface.
%     config_str: string defining the FOV configuration, used by get_fov_param.m
%                 or a config struct compatible with 'get_fov_param.m'
%                 output.
%           camera_fov_param.fov_diag : diagonal field of view angle;
%           camera_fov_param.fov_x : x (horizontal) field of view angle;
%           camera_fov_param.fov_y : y (vertical) field of view angle;
%           camera_fov_param.size_X = frame pixel size in x (horizontal);
%           camera_fov_param.size_Y = frame pixel size in y (vertical);
%           camera_fov_param.camera_offset_Z = vertical (Z) distance in meters; 
%
%   Output:
%     X: 2D matrix (Ny,Nx) containing the distance in meters to camera 
%        center in horizontal direction (x direction). Grid monotonic
%        ascending.
%     Y: 2D matrix (Ny,Nx) containing the distance in meters to camera 
%        center in vertical direction (Y direction). Grid monotonic
%        descending.
%     IMG: 2D Intensity image (Ny,Nx). 
%
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 90 deg downwards %
%
%Sony  EXMOR 1/2.3”
%
% %size: 7.664 mm (H) x 5.984 mm (V)
% Professional Lens
%
% The Phantom 3 Professional’s camera has an f/2.8 lens with a 94⁰ field of view,
% virtually eliminating the unwanted lens distortion common in other
% cameras not built to shoot from the air. Built from 9 separate elements,
% including two aspherical elements, that reduce weight and complexity
% without sacrificing image quality, this lens brings the brightest,
% truest colors to life.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get camera configuration
size_X = size(img,2);
size_Y = size(img,1);

if isstring(config_str) % check if is a string defining the configuration
    camera_fov_param = get_fov_param(config_str,size_X,size_Y);
elseif isstruct(config_str) && isfield(config_str,'size_X') && isfield(config_str,'size_Y')
    
    camera_fov_param = config_str;
    
    % check correct size_X
    if camera_fov_param.size_X ~= size_X 
        error('gridDJIFrame_nadir: config_str wrong size X');
    end
    % check correct size_Y
    if camera_fov_param.size_Y ~= size_Y 
        error('gridDJIFrame_nadir: config_str wrong size Y');
    end
        
else
    error('gridDJIFrame_nadir: config_str config unknown');
end

altitude_after_offset =  (altitude + camera_fov_param.camera_offset_Z);
% altitude_after_offset = -(altitude + camera_fov_param.camera_offset_Z);

x_vec =  linspace(altitude_after_offset*tand(-camera_fov_param.fov_x/2),altitude_after_offset*tand(camera_fov_param.fov_x/2),size_X);
y_vec =  linspace(altitude_after_offset*tand( camera_fov_param.fov_y/2),altitude_after_offset*tand(-camera_fov_param.fov_y/2),size_Y);
% y_vec =  linspace(altitude_after_offset*tand(-camera_fov_param.fov_y/2),altitude_after_offset*tand(camera_fov_param.fov_y/2),size_Y);

[X,Y] = meshgrid(x_vec,y_vec);
IMG = img;

end
